#!/bin/bash
#SBATCH --job-name=fst_chr
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=16:00:00
#SBATCH --output=logs/fst_chr_%A_%a.log

# =============================================================================
# p06e_fst_bychrom.sh — per-chromosome pairwise FST, avoids the OOM that hit
# genome-wide 2D-SFS estimation even at 200G (the same memory-scales-with-
# total-sites lesson as p06d_theta_bychrom.sh, extended to two populations
# estimated jointly). Reuses the existing per-chromosome SAF parts in
# saf_parts/ — no new SAF generation needed, they're already the
# BAQ/C/uniqueOnly-corrected ones.
#
# Modes:
#   --sfs2d-array  (ARRAY TASK) 2D-SFS on one chromosome's SAF pair
#   --sfs2d-sum    sum per-chromosome 2D-SFS vectors -> genome-wide 2D-SFS
#   (default)      (ARRAY TASK) fst index + windowed stats on one chromosome,
#                  using the genome-wide 2D-SFS as the prior
#   --merge        concatenate per-chromosome windowed FST + compute a
#                  sites-weighted global FST summary
#
# USAGE (per pair P1 P2):
#   NCHR=$(wc -l < chromosomes.txt)
#   sbatch --array=1-$NCHR%13 --export=ALL,POP1=pop_map_CAI,POP2=pop_map_NS,OUT_DIR=... \
#     p06e_fst_bychrom.sh --sfs2d-array
#   sbatch --wrap="bash p06e_fst_bychrom.sh --sfs2d-sum POP1=pop_map_CAI POP2=pop_map_NS OUT_DIR=..."
#   sbatch --array=1-$NCHR%13 --export=ALL,POP1=pop_map_CAI,POP2=pop_map_NS,OUT_DIR=... \
#     p06e_fst_bychrom.sh
#   sbatch --wrap="bash p06e_fst_bychrom.sh --merge POP1=pop_map_CAI POP2=pop_map_NS OUT_DIR=..."
# =============================================================================

# ---- mode: --sfs2d-array (ARRAY TASK: one chromosome, joint 2D-SFS) --------
if [[ "${1:-}" == "--sfs2d-array" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    : "${POP1:?POP1 required}" ; : "${POP2:?POP2 required}" ; : "${OUT_DIR:?OUT_DIR required}"
    source ~/.bash_profile
    conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
    set -eu
    THREADS=${SLURM_CPUS_PER_TASK:-4}
    CHROM_FILE=${CHROM_FILE:-$OUT_DIR/chromosomes.txt}
    PAIR="${POP1#pop_map_}_${POP2#pop_map_}"

    CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_FILE")
    [[ -n "$CHR" ]] || { echo "ERROR: no chromosome on line $SLURM_ARRAY_TASK_ID"; exit 1; }
    IDX1="$OUT_DIR/saf_parts/${POP1}.${CHR}.saf.idx"
    IDX2="$OUT_DIR/saf_parts/${POP2}.${CHR}.saf.idx"
    [[ -f "$IDX1" ]] || { echo "ERROR: missing $IDX1"; exit 1; }
    [[ -f "$IDX2" ]] || { echo "ERROR: missing $IDX2"; exit 1; }

    mkdir -p "$OUT_DIR/fst_corrected"
    OUT="$OUT_DIR/fst_corrected/${PAIR}.${CHR}.2dsfs"
    [[ -s "$OUT" ]] && { echo "Already done: $OUT"; exit 0; }

    echo "2D-SFS: $PAIR $CHR (threads=$THREADS)"
    realSFS "$IDX1" "$IDX2" -fold 1 -P "$THREADS" -nSites 2000000 > "$OUT" 2> "${OUT}.log"
    [[ -s "$OUT" ]] || { echo "ERROR: 2D-SFS failed for $CHR"; cat "${OUT}.log"; exit 1; }
    echo "Done: $OUT"
    exit 0
fi

# ---- mode: --sfs2d-sum (sum per-chromosome 2D-SFS -> genome-wide) ---------
if [[ "${1:-}" == "--sfs2d-sum" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    : "${POP1:?POP1 required}" ; : "${POP2:?POP2 required}" ; : "${OUT_DIR:?OUT_DIR required}"
    set -eu
    PAIR="${POP1#pop_map_}_${POP2#pop_map_}"

    PARTS=( "$OUT_DIR/fst_corrected/${PAIR}."*.2dsfs )
    [[ -e "${PARTS[0]}" ]] || { echo "ERROR: no per-chromosome 2D-SFS found for $PAIR"; exit 1; }
    echo "Summing ${#PARTS[@]} per-chromosome 2D-SFS vectors for $PAIR..."

    # realSFS with -nSites prints ONE ROW PER EM CHUNK, not one row per file —
    # collapse each chromosome's rows into a single per-chromosome vector first
    # (element-wise sum within the file), THEN sum those per-chromosome vectors
    # across all 13 chromosomes. Doing this in one flat pass (summing every row
    # of every file together) silently mixes chunk-rows from different files as
    # if they were the same reduction level and produces meaningless output.
    TMP_DIR=$(mktemp -d)
    for f in "${PARTS[@]}"; do
        cname=$(basename "$f")
        awk '{
            n=NF
            for(i=1;i<=n;i++) sum[i]+=$i
        } END {
            out=""
            for(i=1;i<=n;i++) out = out (i>1?" ":"") sum[i]
            print out
        }' "$f" > "$TMP_DIR/$cname"
    done
    awk '{
        n=NF
        for(i=1;i<=n;i++) sum[i]+=$i
    } END {
        out=""
        for(i=1;i<=n;i++) out = out (i>1?" ":"") sum[i]
        print out
    }' "$TMP_DIR"/* > "$OUT_DIR/fst_corrected/${PAIR}.2dsfs.tmp"
    rm -rf "$TMP_DIR"
    mv "$OUT_DIR/fst_corrected/${PAIR}.2dsfs.tmp" "$OUT_DIR/fst_corrected/${PAIR}.2dsfs"

    [[ -s "$OUT_DIR/fst_corrected/${PAIR}.2dsfs" ]] || { echo "ERROR: summed 2D-SFS empty"; exit 1; }
    echo "Genome-wide 2D-SFS written -> $OUT_DIR/fst_corrected/${PAIR}.2dsfs"
    exit 0
fi

# ---- mode: --merge (concatenate windowed FST + weighted global summary) ---
if [[ "${1:-}" == "--merge" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    : "${POP1:?POP1 required}" ; : "${POP2:?POP2 required}" ; : "${OUT_DIR:?OUT_DIR required}"
    set -eu
    PAIR="${POP1#pop_map_}_${POP2#pop_map_}"

    PARTS=( "$OUT_DIR/fst_corrected/${PAIR}."*.fst.windowed.tsv )
    [[ -e "${PARTS[0]}" ]] || { echo "ERROR: no per-chromosome windowed FST found for $PAIR"; exit 1; }

    OUT="$OUT_DIR/fst_corrected/${PAIR}.fst.windowed.tsv"
    head -1 "${PARTS[0]}" > "$OUT"
    for f in "${PARTS[@]}"; do tail -n +2 "$f" >> "$OUT"; done
    N=$(( $(wc -l < "$OUT") - 1 ))
    echo "Merged ${#PARTS[@]} chromosomes -> $OUT ($N windows)"

    # sites-weighted global FST from the windowed output (columns: region chr midPos Nsites Fst)
    awk -F'\t' 'NR>1 && $4>0 {sumFN += $5*$4; sumN += $4} END{
        if (sumN>0) printf "sites_weighted_genomewide_Fst\t%.6f\ttotal_sites\t%d\n", sumFN/sumN, sumN
    }' "$OUT" > "$OUT_DIR/fst_corrected/${PAIR}.fst.global.txt"

    echo "Global (sites-weighted) summary:"
    cat "$OUT_DIR/fst_corrected/${PAIR}.fst.global.txt"
    exit 0
fi

# ---- default mode: per-chromosome array task (fst index + windowed stats) -
_SAVED_ARGS=("$@")
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -- "${_SAVED_ARGS[@]}"
set -eu

: "${POP1:?POP1 required}" ; : "${POP2:?POP2 required}"
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
CHROM_FILE=${CHROM_FILE:-$OUT_DIR/chromosomes.txt}
WINDOW_SIZE=${WINDOW_SIZE:-50000}
WINDOW_STEP=${WINDOW_STEP:-25000}
THREADS=${SLURM_CPUS_PER_TASK:-4}
PAIR="${POP1#pop_map_}_${POP2#pop_map_}"

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_FILE")
[[ -n "$CHR" ]] || { echo "ERROR: no chromosome on line $SLURM_ARRAY_TASK_ID"; exit 1; }

IDX1="$OUT_DIR/saf_parts/${POP1}.${CHR}.saf.idx"
IDX2="$OUT_DIR/saf_parts/${POP2}.${CHR}.saf.idx"
GENOME_2DSFS="$OUT_DIR/fst_corrected/${PAIR}.2dsfs"
[[ -f "$IDX1" ]]          || { echo "ERROR: missing $IDX1"; exit 1; }
[[ -f "$IDX2" ]]          || { echo "ERROR: missing $IDX2"; exit 1; }
[[ -s "$GENOME_2DSFS" ]]  || { echo "ERROR: missing genome-wide 2D-SFS $GENOME_2DSFS (run --sfs2d-sum first)"; exit 1; }

OUT_PREFIX="$OUT_DIR/fst_corrected/${PAIR}.${CHR}"
if [[ -s "${OUT_PREFIX}.fst.windowed.tsv" ]]; then
    echo "Already done: ${OUT_PREFIX}.fst.windowed.tsv"; exit 0
fi

echo "FST  pair=$PAIR  chr=$CHR  threads=$THREADS"

realSFS fst index "$IDX1" "$IDX2" -sfs "$GENOME_2DSFS" -fold 1 -fstout "$OUT_PREFIX" -whichFst 1
[[ -f "${OUT_PREFIX}.fst.idx" ]] || { echo "ERROR: fst.idx not created for $CHR"; exit 1; }

realSFS fst stats2 "${OUT_PREFIX}.fst.idx" -win "$WINDOW_SIZE" -step "$WINDOW_STEP" > "${OUT_PREFIX}.fst.windowed.tsv"
[[ -s "${OUT_PREFIX}.fst.windowed.tsv" ]] || { echo "ERROR: windowed FST not created for $CHR"; exit 1; }

echo "Done: ${OUT_PREFIX}.fst.windowed.tsv"
