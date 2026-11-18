#!/bin/bash
#SBATCH --job-name=beagle_chr
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/beagle_chr_%A_%a.log

# =============================================================================
# p06c2_beagle_bychrom.sh — generate the ALL-samples beagle for PCAngsd,
# one chromosome per array task, then merge and run PCAngsd once.
#
# WHY: the genome-wide `angsd -doGlf 2 ... -sites LDpruned_snps.list` call
# still has to READ every base of every chromosome to check site membership,
# even though it only EMITS ~13k sites. That full-genome node pool at -P 4-40
# is what OOMs (observed at both 40 and 8 threads). Splitting by chromosome
# caps memory to one chromosome and lets all 13 run in parallel.
#
# USAGE:
#   NCHR=$(wc -l < <dir>/chromosomes.txt)
#   sbatch --array=1-$NCHR --export=ALL,BAM_LIST=<list>,REFERENCE=<ref>,OUT_DIR=<dir> \
#     scripts/p06c2_beagle_bychrom.sh
#   # once the array COMPLETES:
#   bash scripts/p06c2_beagle_bychrom.sh --merge OUT_DIR=<dir> PREFIX=pcangsd_ALL
#   # then run PCAngsd on the merged beagle (see --pcangsd mode)
#   sbatch scripts/p06c2_beagle_bychrom.sh --pcangsd OUT_DIR=<dir> PREFIX=pcangsd_ALL THREADS=8
# =============================================================================

# ---- merge mode --------------------------------------------------------------
if [[ "${1:-}" == "--merge" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    : "${OUT_DIR:?OUT_DIR required}"; PREFIX=${PREFIX:-pcangsd_ALL}
    set -eu

    PARTS=( "$OUT_DIR/beagle_parts/${PREFIX}."*.beagle.gz )
    [[ -e "${PARTS[0]}" ]] || { echo "ERROR: no beagle parts found in $OUT_DIR/beagle_parts"; exit 1; }
    echo "Merging ${#PARTS[@]} per-chromosome beagle files..."

    OUT="$OUT_DIR/${PREFIX}.beagle.gz"
    # header from the first part, then all data rows (no header) from every part
    { zcat "${PARTS[0]}" | head -1
      for f in "${PARTS[@]}"; do zcat "$f" | tail -n +2; done
    } | gzip > "$OUT"

    N=$(( $(zcat "$OUT" | wc -l) - 1 ))
    echo "Merged -> $OUT ($N sites)"
    exit 0
fi

# ---- pcangsd mode: run PCAngsd on the already-merged beagle -----------------
if [[ "${1:-}" == "--pcangsd" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    : "${OUT_DIR:?OUT_DIR required}"; PREFIX=${PREFIX:-pcangsd_ALL}
    source ~/.bash_profile
    conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
    set -eu
    THREADS=${THREADS:-8}
    PCANGSD_BIN=${PCANGSD_BIN:-/mnt/parscratch/users/bi4og/conda_envs/pcangsd310/bin/pcangsd}

    BEAGLE="$OUT_DIR/${PREFIX}.beagle.gz"
    [[ -s "$BEAGLE" ]] || { echo "ERROR: merged beagle not found: $BEAGLE"; exit 1; }

    echo "Running PCAngsd on $BEAGLE ..."
    "$PCANGSD_BIN" --beagle "$BEAGLE" -o "$OUT_DIR/$PREFIX" --threads "$THREADS"
    [[ -s "$OUT_DIR/$PREFIX.cov" ]] || { echo "ERROR: PCAngsd did not produce .cov"; exit 1; }
    echo "Done: $OUT_DIR/$PREFIX.cov"
    exit 0
fi

# ---- default mode: per-chromosome array task (ANGSD beagle) ----------------
_SAVED_ARGS=("$@")
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -- "${_SAVED_ARGS[@]}"
set -eu

: "${BAM_LIST:?BAM_LIST required}" ; : "${REFERENCE:?REFERENCE required}"
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
PREFIX=${PREFIX:-pcangsd_ALL}
CHROM_FILE=${CHROM_FILE:-$OUT_DIR/chromosomes.txt}
SITES=${SITES:-$OUT_DIR/LDpruned_snps.list}
MIN_MAPQ=${MIN_MAPQ:-30}
MIN_Q=${MIN_Q:-20}
THREADS=${SLURM_CPUS_PER_TASK:-2}

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_FILE")
[[ -n "$CHR" ]] || { echo "ERROR: no chromosome on line $SLURM_ARRAY_TASK_ID"; exit 1; }

mkdir -p "$OUT_DIR/beagle_parts"
OUT_PREFIX="$OUT_DIR/beagle_parts/${PREFIX}.${CHR}"

if [[ -s "$OUT_PREFIX.beagle.gz" ]]; then
    echo "Already done: $OUT_PREFIX.beagle.gz"; exit 0
fi

echo "=================================================================="
echo "Beagle  chr=$CHR  threads=$THREADS"
echo "=================================================================="

angsd -bam "$BAM_LIST" -ref "$REFERENCE" \
    -r "$CHR" \
    -GL 2 -doGlf 2 -doMajorMinor 3 \
    -sites "$SITES" \
    -doCounts 1 -doMaf 1 \
    -minMapQ "$MIN_MAPQ" -minQ "$MIN_Q" \
    -nQueueSize 50 -P "$THREADS" \
    -out "$OUT_PREFIX"

[[ -s "$OUT_PREFIX.beagle.gz" ]] || { echo "ERROR: beagle not created for $CHR"; exit 1; }
echo "Done: $OUT_PREFIX.beagle.gz"
