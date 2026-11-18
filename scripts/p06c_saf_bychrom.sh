#!/bin/bash
#SBATCH --job-name=saf_chr
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/saf_chr_%A_%a.log

# =============================================================================
# p06c_saf_bychrom.sh — ALL-SITES SAF generation, one chromosome per array task.
#
# WHY: a single genome-wide ANGSD SAF run holds the whole node pool in memory
# and OOMs, which is why -sites subsampling (~0.7% of sites) was used before.
# Restricting each task to one chromosome (-r) caps memory at that chromosome's
# footprint, so every site can be used. The per-chromosome SAFs are then merged
# losslessly with `realSFS cat`.
#
# USAGE — submit one array per population:
#   NCHR=$(wc -l < physalia/angsd_results/chromosomes.txt)
#   sbatch --array=1-$NCHR --job-name=saf_SI \
#     --export=ALL,POP=pop_map_SI,BAM_LIST=<list>,REFERENCE=<ref>,OUT_DIR=<dir> \
#     scripts/p06c_saf_bychrom.sh
#
# then, once the array COMPLETES:
#   bash scripts/p06c_saf_bychrom.sh --merge POP=pop_map_SI OUT_DIR=<dir>
# =============================================================================

# ---- merge mode (run on the login node; it is cheap) ------------------------
if [[ "${1:-}" == "--merge" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    : "${POP:?POP required}" ; : "${OUT_DIR:?OUT_DIR required}"
    source ~/.bash_profile
    conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
    set -eu

    SAF_DIR=${SAF_DIR:-saf_parts}
    PARTS=( "$OUT_DIR/$SAF_DIR/${POP}."*.saf.idx )
    [[ -e "${PARTS[0]}" ]] || { echo "ERROR: no parts found in $OUT_DIR/$SAF_DIR for $POP"; exit 1; }
    echo "Merging ${#PARTS[@]} per-chromosome SAFs for $POP"
    realSFS cat "${PARTS[@]}" -outnames "$OUT_DIR/${POP}"
    [[ -s "$OUT_DIR/${POP}.saf.idx" ]] || { echo "ERROR: merge produced no index"; exit 1; }
    echo "Merged -> $OUT_DIR/${POP}.saf.idx"
    echo "First contig: $(realSFS print "$OUT_DIR/${POP}.saf.idx" 2>/dev/null | head -1 | cut -f1)"
    exit 0
fi

# ---- per-chromosome mode (array task) ---------------------------------------
_SAVED_ARGS=("$@")
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -- "${_SAVED_ARGS[@]}"
set -eu

: "${POP:?POP required (e.g. pop_map_SI)}"
: "${BAM_LIST:?BAM_LIST required}"
: "${REFERENCE:?REFERENCE required}"
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
CHROM_FILE=${CHROM_FILE:-$OUT_DIR/chromosomes.txt}
MIN_MAPQ=${MIN_MAPQ:-20}
MIN_Q=${MIN_Q:-20}
MIN_IND_RATIO=${MIN_IND_RATIO:-0.5}
MAX_DEPTH_IND=${MAX_DEPTH_IND:-50}
THREADS=${SLURM_CPUS_PER_TASK:-4}

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_FILE")
[[ -n "$CHR" ]] || { echo "ERROR: no chromosome on line $SLURM_ARRAY_TASK_ID of $CHROM_FILE"; exit 1; }

N_SAMPLES=$(grep -c . "$BAM_LIST")
MIN_IND=$(awk -v n="$N_SAMPLES" -v r="$MIN_IND_RATIO" 'BEGIN{v=int(n*r+0.5); if(v<1)v=1; print v}')

# n<2 guard: SAF/SFS are meaningless for a single sample
if (( N_SAMPLES < 2 )); then
    echo "SKIP: $POP has n=$N_SAMPLES (<2)"; exit 0
fi

mkdir -p "$OUT_DIR/saf_parts"
PREFIX="$OUT_DIR/saf_parts/${POP}.${CHR}"

if [[ -s "$PREFIX.saf.idx" ]]; then
    echo "Already done: $PREFIX.saf.idx"; exit 0
fi

echo "=================================================================="
echo "SAF  pop=$POP  chr=$CHR  n=$N_SAMPLES  minInd=$MIN_IND  threads=$THREADS"
echo "ALL sites (no -sites subsampling)"
echo "=================================================================="

angsd -bam "$BAM_LIST" -ref "$REFERENCE" -anc "$REFERENCE" \
    -r "$CHR" \
    -GL 2 -doSaf 1 -doCounts 1 \
    -minMapQ "$MIN_MAPQ" -minQ "$MIN_Q" -minInd "$MIN_IND" \
    -setMinDepthInd 1 -setMaxDepthInd "$MAX_DEPTH_IND" \
    -baq 1 -C 50 -uniqueOnly 1 \
    -remove_bads 1 -trim 0 -nQueueSize 50 -P "$THREADS" \
    ${EXTRA_FLAGS:-} \
    -out "$PREFIX"

[[ -s "$PREFIX.saf.idx" ]] || { echo "ERROR: SAF not created for $CHR"; exit 1; }
echo "Done: $PREFIX.saf.idx"
