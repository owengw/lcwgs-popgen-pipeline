#!/bin/bash
#SBATCH --job-name=theta_chr
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/theta_chr_%A_%a.log

# =============================================================================
# p06d_theta_bychrom.sh — ALL-SITES theta/Tajima's D, one chromosome per array
# task, using the already-generated per-chromosome SAF parts in saf_parts/.
#
# WHY: p06d_merge_theta.sh runs `realSFS` + `saf2theta` serially over the
# WHOLE genome SAF (tens of GB), which is both memory-heavy AND slow
# (8-9h wall-clock for a single population). Splitting by chromosome:
#   - caps memory to one chromosome's SAF (a few GB, not 50-85GB)
#   - runs all 13 chromosomes in parallel instead of serially
#   - typically finishes in well under an hour of wall-clock per population
#
# The folded/unfolded SFS is a sum of independent per-site likelihood
# contributions, so summing per-chromosome SFS vectors gives the same
# genome-wide SFS as running on the concatenated data (this is the same
# approach used for the SAF `realSFS cat` merge). Per-chromosome thetas are
# computed against that shared genome-wide SFS as the prior, then the
# per-window .pestPG rows are simply concatenated (each chromosome only
# contributes its own windows).
#
# USAGE:
#   1) Compute the genome-wide SFS from parts (once per population, cheap):
#        bash scripts/p06d_theta_bychrom.sh --sfs POP=pop_map_SI OUT_DIR=<dir>
#   2) Submit one array task per chromosome for thetas:
#        NCHR=$(wc -l < <dir>/chromosomes.txt)
#        sbatch --array=1-$NCHR --export=ALL,POP=pop_map_SI,OUT_DIR=<dir> \
#          scripts/p06d_theta_bychrom.sh
#   3) Merge the per-chromosome .pestPG into one genome-wide file:
#        bash scripts/p06d_theta_bychrom.sh --merge POP=pop_map_SI OUT_DIR=<dir>
# =============================================================================

# ---- mode: --sfs-array (ARRAY TASK: one chromosome, SFS on that part only) --
if [[ "${1:-}" == "--sfs-array" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    : "${POP:?POP required}" ; : "${OUT_DIR:?OUT_DIR required}"
    source ~/.bash_profile
    conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
    set -eu
    THREADS=${SLURM_CPUS_PER_TASK:-2}
    FOLD_SFS=${FOLD_SFS:-1}
    CHROM_FILE=${CHROM_FILE:-$OUT_DIR/chromosomes.txt}

    CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_FILE")
    [[ -n "$CHR" ]] || { echo "ERROR: no chromosome on line $SLURM_ARRAY_TASK_ID"; exit 1; }
    SAF_DIR=${SAF_DIR:-saf_parts}
    IDX="$OUT_DIR/$SAF_DIR/${POP}.${CHR}.saf.idx"
    [[ -f "$IDX" ]] || { echo "ERROR: missing $IDX"; exit 1; }

    [[ -s "$OUT_DIR/${POP}.${CHR}.sfs" ]] && { echo "Already done: $CHR"; exit 0; }

    echo "SFS: $POP $CHR (threads=$THREADS)"
    realSFS "$IDX" -fold "$FOLD_SFS" -P "$THREADS" \
        > "$OUT_DIR/${POP}.${CHR}.sfs" 2>"$OUT_DIR/${POP}.${CHR}.sfs.log"
    [[ -s "$OUT_DIR/${POP}.${CHR}.sfs" ]] || { echo "ERROR: SFS failed for $CHR"; exit 1; }
    echo "Done: $OUT_DIR/${POP}.${CHR}.sfs"
    exit 0
fi

# ---- mode: --sfs-sum (cheap, run AFTER the --sfs-array array completes) ----
if [[ "${1:-}" == "--sfs-sum" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    : "${POP:?POP required}" ; : "${OUT_DIR:?OUT_DIR required}"
    set -eu

    SFS_PARTS=( "$OUT_DIR/${POP}."*.sfs )
    [[ -e "${SFS_PARTS[0]}" ]] || { echo "ERROR: no per-chromosome SFS found for $POP"; exit 1; }
    echo "Summing ${#SFS_PARTS[@]} per-chromosome SFS vectors for $POP..."

    # sum the per-chromosome SFS vectors element-wise -> genome-wide SFS
    awk '{
        n=NF
        for(i=1;i<=n;i++) sum[i]+=$i
    } END {
        out=""
        for(i=1;i<=n;i++) out = out (i>1?" ":"") sum[i]
        print out
    }' "${SFS_PARTS[@]}" > "$OUT_DIR/${POP}.sfs.tmp"
    mv "$OUT_DIR/${POP}.sfs.tmp" "$OUT_DIR/${POP}.sfs"

    [[ -s "$OUT_DIR/${POP}.sfs" ]] || { echo "ERROR: summed genome-wide SFS is empty"; exit 1; }
    echo "Genome-wide SFS written -> $OUT_DIR/${POP}.sfs"
    exit 0
fi

# ---- mode: --merge (concatenate per-chromosome pestPG into genome-wide) ----
if [[ "${1:-}" == "--merge" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    : "${POP:?POP required}" ; : "${OUT_DIR:?OUT_DIR required}"
    set -eu

    PARTS=( "$OUT_DIR/theta_parts/${POP}."*.pestPG )
    [[ -e "${PARTS[0]}" ]] || { echo "ERROR: no per-chromosome pestPG found for $POP"; exit 1; }

    OUT="$OUT_DIR/${POP}.pestPG"
    head -1 "${PARTS[0]}" > "$OUT"
    for f in "${PARTS[@]}"; do tail -n +2 "$f" >> "$OUT"; done
    N=$(($(wc -l < "$OUT") - 1))
    echo "Merged ${#PARTS[@]} chromosomes -> $OUT ($N windows total)"
    exit 0
fi

# ---- default mode: per-chromosome array task (saf2theta) -------------------
_SAVED_ARGS=("$@")
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -- "${_SAVED_ARGS[@]}"
set -eu

: "${POP:?POP required}"
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
CHROM_FILE=${CHROM_FILE:-$OUT_DIR/chromosomes.txt}
WINDOW_SIZE=${WINDOW_SIZE:-50000}
WINDOW_STEP=${WINDOW_STEP:-10000}
THREADS=${SLURM_CPUS_PER_TASK:-2}

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_FILE")
[[ -n "$CHR" ]] || { echo "ERROR: no chromosome on line $SLURM_ARRAY_TASK_ID"; exit 1; }

SAF_DIR=${SAF_DIR:-saf_parts}
SAF_IDX="$OUT_DIR/$SAF_DIR/${POP}.${CHR}.saf.idx"
GENOME_SFS="$OUT_DIR/${POP}.sfs"
[[ -f "$SAF_IDX" ]]     || { echo "ERROR: missing $SAF_IDX"; exit 1; }
[[ -s "$GENOME_SFS" ]]  || { echo "ERROR: missing genome-wide SFS $GENOME_SFS (run --sfs first)"; exit 1; }

mkdir -p "$OUT_DIR/theta_parts"
OUT_PREFIX="$OUT_DIR/theta_parts/${POP}.${CHR}"

if [[ -s "$OUT_PREFIX.pestPG" ]]; then
    echo "Already done: $OUT_PREFIX.pestPG"; exit 0
fi

echo "=================================================================="
echo "Theta  pop=$POP  chr=$CHR  threads=$THREADS"
echo "=================================================================="

realSFS saf2theta "$SAF_IDX" -sfs "$GENOME_SFS" -outname "$OUT_PREFIX" -P "$THREADS" -fold 1
[[ -f "$OUT_PREFIX.thetas.idx" ]] || { echo "ERROR: thetas.idx not created for $CHR"; exit 1; }

thetaStat do_stat "$OUT_PREFIX.thetas.idx" -win "$WINDOW_SIZE" -step "$WINDOW_STEP" -outnames "$OUT_PREFIX"
[[ -s "$OUT_PREFIX.pestPG" ]] || { echo "ERROR: pestPG not created for $CHR"; exit 1; }

echo "Done: $OUT_PREFIX.pestPG"
