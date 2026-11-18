#!/bin/bash
#SBATCH --job-name=ind_het
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=1-166%40
#SBATCH --output=logs/p07a_ind_het_%A_%a.log

# =============================================================================
# p07a_individual_het.sh — RagTag rewrite, ALL SITES (no random-10M subsample).
# Per-individual 1-sample folded SFS heterozygosity. One array task per
# individual (180, not 198/183). No EXCLUDE regex — reads directly from the
# final 180-sample bam_list_dedup_overlapclipped.list, so every excluded/
# replicate sample from the session's QC is already correctly absent.
#
# Memory note: a single-sample genome-wide SAF is far lighter than the
# population-level ones that needed per-chromosome splitting (n=1 vs n=16-85),
# so one job per individual, restricted to the 13 placed chromosomes via
# -rf, is safe at 16G.
# =============================================================================

_SAVED_ARGS=("$@")
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -- "${_SAVED_ARGS[@]}"
set -eu

OUT_DIR="physalia/angsd_results"
REFERENCE="/mnt/parscratch/users/bi4og/ncbi_submission/plongirostris.fasta"
CHROM_FILE="$OUT_DIR/chromosomes.txt"
BAM_LIST="physalia/sample_lists/bam_list_dedup_overlapclipped.list"
HET_DIR="$OUT_DIR/heterozygosity_corrected"
THREADS=${SLURM_CPUS_PER_TASK:-2}

mkdir -p "$HET_DIR" logs

TOTAL=$(wc -l < "$BAM_LIST")
if [[ ${SLURM_ARRAY_TASK_ID} -gt $TOTAL ]]; then
    echo "Array task ${SLURM_ARRAY_TASK_ID} > $TOTAL samples — exiting"
    exit 0
fi

BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_LIST")
[[ -n "$BAM" ]] || { echo "ERROR: no BAM on line $SLURM_ARRAY_TASK_ID"; exit 1; }
SAMPLE_ID=$(basename "$BAM" | sed 's/_merged.*//')
echo "Processing: $SAMPLE_ID (task $SLURM_ARRAY_TASK_ID / $TOTAL)"

SAF_PREFIX="$HET_DIR/${SAMPLE_ID}"
SFS_FILE="${SAF_PREFIX}.sfs"
HET_FILE="${SAF_PREFIX}.het"

if [[ -s "$HET_FILE" ]]; then
    echo "Already done: $HET_FILE"
    cat "$HET_FILE"
    exit 0
fi

echo "=================================================================="
echo "ANGSD SAF (all sites, all 13 placed chromosomes)"
echo "=================================================================="
angsd \
    -i "$BAM" \
    -ref "$REFERENCE" \
    -anc "$REFERENCE" \
    -rf "$CHROM_FILE" \
    -GL 2 \
    -doSaf 1 \
    -doCounts 1 \
    -minMapQ 20 \
    -minQ 20 \
    -baq 1 \
    -C 50 \
    -uniqueOnly 1 \
    -setMinDepthInd 1 \
    -setMaxDepthInd 50 \
    -remove_bads 1 \
    -trim 0 \
    -P $THREADS \
    -nQueueSize 50 \
    -out "$SAF_PREFIX" \
    2> "${SAF_PREFIX}.angsd.log"

[[ -s "${SAF_PREFIX}.saf.idx" ]] || { echo "ERROR: SAF not created for $SAMPLE_ID"; tail -30 "${SAF_PREFIX}.angsd.log"; exit 1; }

realSFS \
    "${SAF_PREFIX}.saf.idx" \
    -fold 1 \
    -P $THREADS \
    > "$SFS_FILE" \
    2> "${SAF_PREFIX}.sfs.log"

[[ -s "$SFS_FILE" ]] || { echo "ERROR: realSFS produced no output for $SAMPLE_ID"; cat "${SAF_PREFIX}.sfs.log"; exit 1; }

# Extract heterozygosity: het / (hom_ref + het + hom_alt)
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

echo "Complete: $SAMPLE_ID — het=$(cut -f2 "$HET_FILE" | tail -1)"
