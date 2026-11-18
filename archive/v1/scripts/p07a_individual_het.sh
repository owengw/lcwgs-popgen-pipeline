#!/bin/bash
#SBATCH --job-name=ind_het
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-198
#SBATCH --output=logs/p07a_ind_het_%A_%a.log

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -euo pipefail
set -x

# ============================================================
# Variables — edit these if needed
# ============================================================
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
BAM_LIST=${BAM_LIST:-physalia/angsd_results/all_samples.list}
REFERENCE="genome/Plong_genome_flye.fasta"
SITES_FILE="$OUT_DIR/random_10M_sites.list"
RF_FILE="$OUT_DIR/random_10M_contigs.txt"
HET_DIR="$OUT_DIR/heterozygosity_corrected"
THREADS=4

mkdir -p "$HET_DIR"

# ============================================================
# Get BAM file for this array task
# ============================================================
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_LIST")

if [[ -z "$BAM" ]]; then
    echo "ERROR: No BAM found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Extract sample ID from BAM path
SAMPLE_ID=$(basename "$BAM" | sed 's/_merged.*//')
echo "Processing: $SAMPLE_ID (task $SLURM_ARRAY_TASK_ID)"

SAF_PREFIX="$HET_DIR/${SAMPLE_ID}"
SFS_FILE="${SAF_PREFIX}.sfs"
HET_FILE="${SAF_PREFIX}.het"

# ============================================================
# Step 1 — Generate per-individual SAF
# Uses same parameters as population SAF runs:
#   -GL 2 (GATK model)
#   -minMapQ 10, -minQ 20
#   -setMinDepthInd 1 (at least 1 read)
#   -remove_bads 1, -trim 0
#   random_10M_sites for unbiased genome-wide estimate
# Note: no -minInd filter for single-individual runs
# ============================================================
angsd \
    -i "$BAM" \
    -ref "$REFERENCE" \
    -anc "$REFERENCE" \
    -GL 2 \
    -doSaf 1 \
    -doCounts 1 \
    -minMapQ 10 \
    -minQ 20 \
    -setMinDepthInd 1 \
    -remove_bads 1 \
    -trim 0 \
    -P $THREADS \
    -nQueueSize 50 \
    -sites "$SITES_FILE" \
    -rf "$RF_FILE" \
    -out "$SAF_PREFIX"

# ============================================================
# Step 2 — Estimate 1-sample folded SFS using realSFS
# The folded SFS is appropriate when using reference as ancestral
# For a diploid individual: SFS has 3 entries [n_hom_ref, n_het, n_hom_alt]
# ============================================================
realSFS \
    "${SAF_PREFIX}.saf.idx" \
    -fold 1 \
    -P $THREADS \
    > "$SFS_FILE"

# ============================================================
# Step 3 — Extract heterozygosity from SFS
# SFS entries: [count_hom_ref, count_het, count_hom_alt]
# Heterozygosity = count_het / (count_hom_ref + count_het + count_hom_alt)
# This is the correct genome-wide per-site heterozygosity
# ============================================================
# Extract heterozygosity from SFS using awk (no Python needed)
awk -v sample="$SAMPLE_ID" -v out="$HET_FILE" '
BEGIN { print "sample_id\theterozygosity\thom_ref_prop\thom_alt_prop\t" \
              "n_hom_ref\tn_het\tn_hom_alt\tn_sites_total" > out }
{
    hom_ref = $1; het = $2; hom_alt = $3
    total = hom_ref + het + hom_alt
    if (total > 0) {
        heterozygosity = het / total
        printf "%s\t%.8f\t%.8f\t%.8f\t%.1f\t%.1f\t%.1f\t%.1f\n",
               sample, heterozygosity, hom_ref/total, hom_alt/total,
               hom_ref, het, hom_alt, total >> out
        printf "%s: heterozygosity = %.6f (het=%.0f, total=%.0f)\n",
               sample, heterozygosity, het, total
    } else {
        print "ERROR: total sites = 0 for " sample > "/dev/stderr"
        exit 1
    }
}' "$SFS_FILE"

echo "Complete: $SAMPLE_ID — het=$(cut -f2 $HET_FILE | tail -1)"
