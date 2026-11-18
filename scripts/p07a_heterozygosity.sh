#!/bin/bash
#SBATCH --job-name=heterozygosity
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --output=logs/p07a_heterozygosity_%j.log

# Load environment
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -eu
set -x

# Variables
BAM_LIST=${BAM_LIST:-all_bams.list}
REFERENCE=${REFERENCE:-}
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
THREADS=${THREADS:-8}

echo "=================================================================="
echo "Individual Heterozygosity Analysis"
echo "=================================================================="
echo "BAM list: $BAM_LIST"
echo "Output directory: $OUT_DIR"
echo "=================================================================="

mkdir -p ${OUT_DIR}/heterozygosity

# Output prefix
OUT_PREFIX="${OUT_DIR}/heterozygosity/het_all_individuals"

# Generate beagle file from random 10M sites if it doesn't exist
BEAGLE_FILE="${OUT_PREFIX}.beagle.gz"

if [[ ! -f "$BEAGLE_FILE" ]]; then
    echo "Generating beagle file from random 10M sites..."
    
    # Use random sites for unbiased heterozygosity estimation
    SITES_FILE="${OUT_DIR}/random_10M_sites.list"
    RF_FILE="${OUT_DIR}/random_10M_contigs.txt"
    
    if [[ ! -f "$SITES_FILE" ]]; then
        echo "ERROR: Random sites file not found: $SITES_FILE"
        exit 1
    fi
    
    N_SAMPLES=$(wc -l < "$BAM_LIST")
    MIN_IND=$(echo "$N_SAMPLES * 0.30" | bc | awk '{print int($1+0.5)}')
    [[ $MIN_IND -lt 1 ]] && MIN_IND=1
    
    echo "Samples: $N_SAMPLES, Min individuals: $MIN_IND"
    
    angsd -bam "$BAM_LIST" \
        -ref "$REFERENCE" \
        -GL 2 \
        -doGlf 2 \
        -doMajorMinor 1 \
        -doMaf 1 \
        -minMapQ 10 \
        -minQ 20 \
        -minInd $MIN_IND \
        -sites "$SITES_FILE" \
        -rf "$RF_FILE" \
        -P $THREADS \
        -out "$OUT_PREFIX" 2>&1 | tee ${OUT_PREFIX}_beagle.log
    
    if [[ ! -f "$BEAGLE_FILE" ]]; then
        echo "ERROR: Failed to generate beagle file"
        exit 1
    fi
    
    echo "Beagle file generated successfully"
else
    echo "Using existing beagle file: $BEAGLE_FILE"
fi

# Calculate per-individual heterozygosity from beagle file
echo "Calculating per-individual heterozygosity..."

Rscript --vanilla - "$BEAGLE_FILE" "$BAM_LIST" "$OUT_PREFIX" <<'REOF'
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
beagle_file <- args[1]
bam_list <- args[2]
out_prefix <- args[3]

cat("Reading beagle file...\n")
beagle <- read_tsv(beagle_file, col_types = cols(), show_col_types = FALSE)

# Extract sample names from header
header_names <- colnames(beagle)
sample_cols <- header_names[-(1:3)]  # Remove marker, allele1, allele2
sample_indices <- seq(1, length(sample_cols), by = 3)
sample_names <- str_replace(sample_cols[sample_indices], "_0$", "")

cat("Found", length(sample_names), "samples\n")

# Read BAM list to get proper sample IDs
bam_samples <- read_tsv(bam_list, col_names = "bam_path", col_types = cols())
bam_samples <- bam_samples %>%
  mutate(sample_id = basename(bam_path) %>% str_replace("\\.bam$", ""))

cat("Processing", nrow(beagle), "sites...\n")

# Calculate heterozygosity for each individual
het_results <- tibble()

for (i in 1:length(sample_names)) {
  if (i %% 20 == 0) cat("  Processing individual", i, "of", length(sample_names), "...\n")
  
  # Get the 3 GL columns for this individual (GL for 0/0, 0/1, 1/1)
  col_idx <- 3 + (i - 1) * 3 + 1
  gl_cols <- beagle[, col_idx:(col_idx + 2)]
  gl_mat <- as.matrix(gl_cols)
  
  # Convert log10-scaled likelihoods to probabilities
  # GL are already probabilities in beagle format (sum to 1 per site)
  # Just use them directly
  probs <- gl_mat
  
  # For sites with valid data (not all 0.333...)
  valid_sites <- rowSums(abs(probs - 0.333333) > 0.01) > 0
  
  if (sum(valid_sites) < 10) {
    cat("  WARNING: Sample", i, "has very few valid sites\n")
    next
  }
  
  # Expected heterozygosity = P(0/1) for each site
  het_per_site <- probs[valid_sites, 2]
  
  # Calculate statistics
  n_sites <- sum(valid_sites)
  mean_het <- mean(het_per_site, na.rm = TRUE)
  median_het <- median(het_per_site, na.rm = TRUE)
  sd_het <- sd(het_per_site, na.rm = TRUE)
  
  # Observed heterozygosity (sites where het prob > 0.5)
  obs_het_sites <- sum(het_per_site > 0.5, na.rm = TRUE)
  obs_het_prop <- obs_het_sites / n_sites
  
  het_results <- bind_rows(het_results, tibble(
    sample_name = sample_names[i],
    sample_id = bam_samples$sample_id[i],
    n_sites = n_sites,
    mean_heterozygosity = mean_het,
    median_heterozygosity = median_het,
    sd_heterozygosity = sd_het,
    observed_het_sites = obs_het_sites,
    observed_het_proportion = obs_het_prop
  ))
}

# Write results
cat("Writing results...\n")
write_tsv(het_results, paste0(out_prefix, "_individual_summary.txt"))

# Summary statistics
cat("\n=== Individual Heterozygosity Summary ===\n")
cat("Number of individuals:", nrow(het_results), "\n")
cat("Mean heterozygosity across all individuals:", mean(het_results$mean_heterozygosity), "\n")
cat("SD:", sd(het_results$mean_heterozygosity), "\n")
cat("Range:", min(het_results$mean_heterozygosity), "-", max(het_results$mean_heterozygosity), "\n")
cat("\nTop 10 most heterozygous individuals:\n")
print(het_results %>% arrange(desc(mean_heterozygosity)) %>% select(sample_id, mean_heterozygosity, n_sites) %>% head(10))
cat("\nBottom 10 least heterozygous individuals:\n")
print(het_results %>% arrange(mean_heterozygosity) %>% select(sample_id, mean_heterozygosity, n_sites) %>% head(10))

# Write overall summary
sink(paste0(out_prefix, "_overall_summary.txt"))
cat("Individual Heterozygosity Summary\n")
cat("==================================\n\n")
cat("Analysis based on random 10M site subsample\n\n")
cat("Number of individuals:", nrow(het_results), "\n")
cat("Sites analyzed per individual: ~", round(mean(het_results$n_sites)), "\n\n")
cat("Mean heterozygosity:\n")
cat("  Across all individuals:", mean(het_results$mean_heterozygosity), "\n")
cat("  Median:", median(het_results$mean_heterozygosity), "\n")
cat("  SD:", sd(het_results$mean_heterozygosity), "\n")
cat("  Range:", min(het_results$mean_heterozygosity), "-", max(het_results$mean_heterozygosity), "\n\n")
cat("Observed heterozygosity proportion:\n")
cat("  Mean:", mean(het_results$observed_het_proportion), "\n")
cat("  Median:", median(het_results$observed_het_proportion), "\n")
cat("  SD:", sd(het_results$observed_het_proportion), "\n")
sink()

cat("\nComplete!\n")
REOF

echo ""
echo "=================================================================="
echo "Heterozygosity Analysis Complete"
echo "=================================================================="
echo "Individual results: ${OUT_PREFIX}_individual_summary.txt"
echo "Overall summary: ${OUT_PREFIX}_overall_summary.txt"
echo "=================================================================="
echo ""
echo "Interpretation:"
echo "  Mean heterozygosity = average probability of being heterozygous per site"
echo "  Values typically range from 0.001 to 0.01 for most species"
echo "  Higher values indicate more genetic diversity within individuals"
echo "=================================================================="