#!/bin/bash
#SBATCH --job-name=rohan_individual
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p08a_rohan_individual_%A_%a.log
#SBATCH --array=1-195

# =============================================================================
# ROHan per-individual ROH estimation
# Runs on each individual BAM file using --tstv mode
#
# --tstv [ratio]: uses a fixed Ts/Tv ratio (default 2.1) rather than requiring
#                 a mutation rate — appropriate for non-model organisms and lcWGS
# NOTE: --size in this ROHan version = window size in bp (default 1,000,000)
#       NOT genome size — genome size is inferred from BAM headers automatically
#       Do not pass genome size as --size
#
# ROHan binary: /mnt/parscratch/users/bi4og/ROHan/bin/rohan
# Requires modules: HTSlib/1.17-GCC-12.2.0 and GSL/2.7-GCC-12.2.0
#
# Array size: 195 (198 total - 3 excluded: SI43, SI45, SI83)
# Verify with: grep -cEv "SI43|SI45|SI83" physalia/angsd_results/all_samples.list
#
# Excludes: SI43, SI45 (failed samples), SI83 (duplicate of SI77)
# o-suffix and T-suffix samples retained — ROHan is per-individual so
# results from duplicates can be verified against each other post-hoc
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
module load HTSlib/1.17-GCC-12.2.0
module load GSL/2.7-GCC-12.2.0

set -euo pipefail

# =============================================================================
# Variables
# =============================================================================
ROHAN="/mnt/parscratch/users/bi4og/ROHan/bin/rohan"
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
REFERENCE="genome/Plong_genome_flye.fasta"
BAM_LIST="$OUT_DIR/all_samples.list"
ROH_DIR="$OUT_DIR/roh/individual"
THREADS=4

EXCLUDE="SI43|SI45|SI83"

mkdir -p "$ROH_DIR"
mkdir -p logs

# =============================================================================
# Build filtered sample list (cached)
# =============================================================================
FILTERED_LIST="$ROH_DIR/rohan_sample.list"
if [[ ! -f "$FILTERED_LIST" ]]; then
    grep -vE "$EXCLUDE" "$BAM_LIST" > "$FILTERED_LIST"
    echo "Filtered sample list created: $(wc -l < "$FILTERED_LIST") samples"
fi

TOTAL=$(wc -l < "$FILTERED_LIST")
echo "Total samples: $TOTAL"

if [[ ${SLURM_ARRAY_TASK_ID} -gt $TOTAL ]]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} exceeds sample count $TOTAL — exiting"
    exit 0
fi

# =============================================================================
# Get this array task's BAM file
# =============================================================================
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILTERED_LIST")
SAMPLE=$(basename "$BAM" | sed 's/_merged.*//')

echo "=================================================================="
echo "ROHan individual: $SAMPLE"
echo "BAM: $BAM"
echo "Array task: ${SLURM_ARRAY_TASK_ID} / $TOTAL"
echo "=================================================================="

if [[ -f "$ROH_DIR/${SAMPLE}.summary" ]]; then
    echo "Output already exists, skipping: $ROH_DIR/${SAMPLE}.summary"
    exit 0
fi

if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM file not found: $BAM"
    exit 1
fi

# =============================================================================
# Run ROHan
# --tstv 2.1: Ts/Tv ratio (2.1 is standard default for vertebrates)
# -t:         number of threads
# -o:         output prefix
# NOTE: --size is window size in bp, not genome size
#       genome size inferred from BAM header automatically
# =============================================================================
"$ROHAN" \
    --tstv 2.1 \
    -t $THREADS \
    -o "$ROH_DIR/${SAMPLE}" \
    "$REFERENCE" \
    "$BAM" \
    2> "$ROH_DIR/${SAMPLE}.log"

if [[ -f "$ROH_DIR/${SAMPLE}.summary" ]]; then
    echo "Done: $SAMPLE"
    cat "$ROH_DIR/${SAMPLE}.summary"
else
    echo "ERROR: No summary file produced — check $ROH_DIR/${SAMPLE}.log"
    tail -10 "$ROH_DIR/${SAMPLE}.log"
    exit 1
fi