#!/bin/bash
#SBATCH --job-name=saf_chr_masked
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/saf_chr_masked_%A_%a.log

# =============================================================================
# p06c_saf_bychrom_masked.sh — REPEAT-EXCLUDED all-sites SAF, one chromosome
# per array task.
#
# WHY: all-sites theta on the unmasked genome gave Tajima's D of -2 to -2.7
# uniformly across every population — implausible as independent biology.
# RepeatMasker (on the same Pilon-polished sequence used here, coordinates
# confirmed identical) shows 52% of the 13 placed chromosomes is repetitive.
# Mismapping in that fraction inflates spurious rare-variant calls, which
# is a far more parsimonious explanation than four populations independently
# showing near-identical strong selection/expansion signals.
#
# This restricts ANGSD to the callable_regions.bed complement (non-repeat)
# via -sites, on top of the existing per-chromosome -r restriction. Same
# memory/time profile as the unmasked version since site count is roughly
# halved, not doubled.
#
# PREREQUISITE (run once, on the login node, cheap):
#   angsd sites index physalia/angsd_results/callable_sites.list
# (this script expects that file/index to already exist — see below)
# =============================================================================

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
CALLABLE_SITES=${CALLABLE_SITES:-$OUT_DIR/callable_sites.list}
MIN_MAPQ=${MIN_MAPQ:-20}
MIN_Q=${MIN_Q:-20}
MIN_IND_RATIO=${MIN_IND_RATIO:-0.5}
MAX_DEPTH_IND=${MAX_DEPTH_IND:-50}
THREADS=${SLURM_CPUS_PER_TASK:-4}

[[ -f "${CALLABLE_SITES}.idx" ]] || { echo "ERROR: $CALLABLE_SITES not indexed. Run: angsd sites index $CALLABLE_SITES"; exit 1; }

CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_FILE")
[[ -n "$CHR" ]] || { echo "ERROR: no chromosome on line $SLURM_ARRAY_TASK_ID"; exit 1; }

N_SAMPLES=$(grep -c . "$BAM_LIST")
MIN_IND=$(awk -v n="$N_SAMPLES" -v r="$MIN_IND_RATIO" 'BEGIN{v=int(n*r+0.5); if(v<1)v=1; print v}')

if (( N_SAMPLES < 2 )); then
    echo "SKIP: $POP has n=$N_SAMPLES (<2)"; exit 0
fi

mkdir -p "$OUT_DIR/saf_parts_masked"
PREFIX="$OUT_DIR/saf_parts_masked/${POP}.${CHR}"

if [[ -s "$PREFIX.saf.idx" ]]; then
    echo "Already done: $PREFIX.saf.idx"; exit 0
fi

echo "=================================================================="
echo "SAF (repeat-masked)  pop=$POP  chr=$CHR  n=$N_SAMPLES  minInd=$MIN_IND"
echo "=================================================================="

angsd -bam "$BAM_LIST" -ref "$REFERENCE" -anc "$REFERENCE" \
    -r "$CHR" \
    -sites "$CALLABLE_SITES" \
    -GL 2 -doSaf 1 -doCounts 1 \
    -minMapQ "$MIN_MAPQ" -minQ "$MIN_Q" -minInd "$MIN_IND" \
    -setMinDepthInd 1 -setMaxDepthInd "$MAX_DEPTH_IND" \
    -remove_bads 1 -trim 0 -nQueueSize 50 -P "$THREADS" \
    -out "$PREFIX"

[[ -s "$PREFIX.saf.idx" ]] || { echo "ERROR: SAF not created for $CHR"; exit 1; }
echo "Done: $PREFIX.saf.idx"
