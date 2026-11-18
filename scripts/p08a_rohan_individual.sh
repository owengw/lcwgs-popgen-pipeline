#!/bin/bash
#SBATCH --job-name=rohan_individual
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p08a_rohan_individual_%A_%a.log
#SBATCH --array=1-166%30

# =============================================================================
# p08a_rohan_individual.sh — RagTag rewrite.
# Per-individual ROH estimation via ROHan, --tstv mode (fixed Ts/Tv ratio,
# appropriate for non-model organisms / lcWGS — no mutation rate needed).
#
# CHANGES FROM ORIGINAL:
#  - array=1-180 (was 1-195), reads bam_list_dedup_overlapclipped.list
#    directly — no EXCLUDE regex. Your final 180-sample list already
#    correctly excludes SI43/SI45/SI35/SB4 and the 14 amplification
#    replicates, while correctly KEEPING SI1T (unique animal, no original).
#  - reference is the real RagTag plongirostris.fasta, not Flye.
#
# NOTE: --size in this ROHan version = window size in bp (default 1Mb),
#       NOT genome size — genome size is inferred from the BAM header.
# ROHan binary: /mnt/parscratch/users/bi4og/ROHan/bin/rohan
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
module load HTSlib/1.17-GCC-12.2.0
module load GSL/2.7-GCC-12.2.0
set -eu

ROHAN="/mnt/parscratch/users/bi4og/ROHan/bin/rohan"
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
REFERENCE="/mnt/parscratch/users/bi4og/ncbi_submission/plongirostris.fasta"
BAM_LIST="physalia/sample_lists/bam_list_dedup_overlapclipped.list"
ROH_DIR="$OUT_DIR/roh/individual"
THREADS=${SLURM_CPUS_PER_TASK:-4}

mkdir -p "$ROH_DIR" logs

TOTAL=$(wc -l < "$BAM_LIST")
if [[ ${SLURM_ARRAY_TASK_ID} -gt $TOTAL ]]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} exceeds sample count $TOTAL — exiting"
    exit 0
fi

BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_LIST")
SAMPLE=$(basename "$BAM" | sed 's/_merged.*//')

echo "=================================================================="
echo "ROHan individual: $SAMPLE"
echo "BAM: $BAM"
echo "Array task: ${SLURM_ARRAY_TASK_ID} / $TOTAL"
echo "=================================================================="

if [[ -s "$ROH_DIR/${SAMPLE}.summary.txt" ]]; then
    echo "Output already exists, skipping: $ROH_DIR/${SAMPLE}.summary.txt"
    exit 0
fi

[[ -f "$BAM" ]] || { echo "ERROR: BAM file not found: $BAM"; exit 1; }

"$ROHAN" \
    --tstv 2.1 \
    -t $THREADS \
    -o "$ROH_DIR/${SAMPLE}" \
    "$REFERENCE" \
    "$BAM" \
    2> "$ROH_DIR/${SAMPLE}.log"

if [[ -s "$ROH_DIR/${SAMPLE}.summary.txt" ]]; then
    echo "Done: $SAMPLE"
    cat "$ROH_DIR/${SAMPLE}.summary.txt"
else
    echo "ERROR: No summary file produced — check $ROH_DIR/${SAMPLE}.log"
    tail -20 "$ROH_DIR/${SAMPLE}.log"
    exit 1
fi
