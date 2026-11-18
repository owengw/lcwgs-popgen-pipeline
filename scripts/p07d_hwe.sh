#!/bin/bash
#SBATCH --job-name=hwe_analysis
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/p07d_hwe_%j.log

# Load environment
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -eu
set -x

# Variables
BAM_LIST=${BAM_LIST:-}
REFERENCE=${REFERENCE:-}
OUT_DIR=${OUT_DIR:-}
THREADS=${THREADS:-16}
MIN_MAF=${MIN_MAF:-0.05}
MIN_DEPTH=${MIN_DEPTH:-10}
MAX_DEPTH=${MAX_DEPTH:-500}

echo "=================================================================="
echo "Hardy-Weinberg Equilibrium Analysis"
echo "=================================================================="
echo "BAM list: $BAM_LIST"
echo "Output directory: $OUT_DIR"
echo "Min MAF: $MIN_MAF"
echo "Min depth per individual: $MIN_DEPTH"
echo "Max depth per individual: $MAX_DEPTH"
echo "=================================================================="

mkdir -p ${OUT_DIR}/hwe

# Extract population name
POP_NAME=$(basename $BAM_LIST .list)
OUT_PREFIX="${OUT_DIR}/hwe/hwe_${POP_NAME}"

# Get sample size
N_SAMPLES=$(wc -l < "$BAM_LIST")
echo "Sample size: $N_SAMPLES"

# Check if population has enough samples (minimum 5 for HWE)
if [[ $N_SAMPLES -lt 5 ]]; then
    echo "WARNING: Population $POP_NAME has only $N_SAMPLES samples"
    echo "HWE analysis requires at least 5 samples, skipping"
    exit 0
fi

# Use random 10M sites for unbiased HWE testing
SITES_FILE="${OUT_DIR}/random_10M_sites.list"
RF_FILE="${OUT_DIR}/random_10M_contigs.txt"

if [[ ! -f "$SITES_FILE" ]]; then
    echo "ERROR: Random sites file not found: $SITES_FILE"
    echo "Run p06c_angsd_subset.sh first to generate random sites"
    exit 1
fi

echo "Using random 10M sites for HWE analysis"

# Run ANGSD HWE analysis
echo "Running ANGSD HWE analysis..."

# Calculate minInd (20% of samples)
MIN_IND=$(echo "$N_SAMPLES * 0.2" | bc | awk '{print int($1+0.5)}')
[[ $MIN_IND -lt 1 ]] && MIN_IND=1

angsd -bam "$BAM_LIST" \
    -ref "$REFERENCE" \
    -GL 1 \
    -doMaf 1 \
    -doMajorMinor 1 \
    -doHWE 2 \
    -minMapQ 10 \
    -minQ 20 \
    -minMaf $MIN_MAF \
    -minInd $MIN_IND \
    -P $THREADS \
    -sites "$SITES_FILE" \
    -rf "$RF_FILE" \
    -out "$OUT_PREFIX" 2>&1 | tee ${OUT_PREFIX}.log

if [[ ! -f "${OUT_PREFIX}.mafs.gz" ]]; then
    echo "ERROR: ANGSD HWE analysis failed"
    tail -50 ${OUT_PREFIX}.log
    exit 1
fi

echo "Extracting F values..."
zcat ${OUT_PREFIX}.mafs.gz | awk 'NR>1 {print $1"\t"$2"\t"$10}' > ${OUT_PREFIX}_F_values.txt

echo "Calculating HWE statistics..."

# Create summary statistics
Rscript --vanilla - "${OUT_PREFIX}_F_values.txt" "${OUT_PREFIX}" <<'REOF'
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
f_file <- args[1]
out_prefix <- args[2]

# Read F values
f_data <- read_tsv(f_file, col_names = c("chr", "pos", "F"), col_types = cols())

cat("Processing", nrow(f_data), "sites with HWE statistics\n")

# Summary statistics
summary_stats <- f_data %>%
  summarise(
    n_sites = n(),
    mean_F = mean(F, na.rm = TRUE),
    median_F = median(F, na.rm = TRUE),
    sd_F = sd(F, na.rm = TRUE),
    min_F = min(F, na.rm = TRUE),
    max_F = max(F, na.rm = TRUE),
    n_positive_F = sum(F > 0, na.rm = TRUE),
    n_negative_F = sum(F < 0, na.rm = TRUE),
    prop_positive_F = sum(F > 0, na.rm = TRUE) / n(),
    n_extreme_positive = sum(F > 0.5, na.rm = TRUE),
    n_extreme_negative = sum(F < -0.5, na.rm = TRUE)
  )

# Write summary
write_tsv(summary_stats, paste0(out_prefix, "_summary.txt"))

# Create histogram data
hist_data <- f_data %>%
  mutate(F_bin = cut(F, breaks = seq(-1, 1, 0.05))) %>%
  count(F_bin) %>%
  drop_na()

write_tsv(hist_data, paste0(out_prefix, "_histogram.txt"))

# Identify extreme sites
extreme_sites <- f_data %>%
  filter(abs(F) > 0.5) %>%
  arrange(desc(abs(F)))

write_tsv(extreme_sites, paste0(out_prefix, "_extreme_F.txt"))

# Print summary
cat("\n=== HWE Summary ===\n")
cat("Total sites:", summary_stats$n_sites, "\n")
cat("Mean F:", summary_stats$mean_F, "\n")
cat("Median F:", summary_stats$median_F, "\n")
cat("SD F:", summary_stats$sd_F, "\n")
cat("Range:", summary_stats$min_F, "to", summary_stats$max_F, "\n")
cat("Positive F (inbreeding):", summary_stats$prop_positive_F * 100, "%\n")
cat("Extreme positive F (>0.5):", summary_stats$n_extreme_positive, "\n")
cat("Extreme negative F (<-0.5):", summary_stats$n_extreme_negative, "\n")
REOF

echo ""
echo "=================================================================="
echo "HWE Analysis Complete"
echo "=================================================================="
echo "Summary: ${OUT_PREFIX}_summary.txt"
echo "F values: ${OUT_PREFIX}_F_values.txt"
echo "Histogram: ${OUT_PREFIX}_histogram.txt"
echo "Extreme sites: ${OUT_PREFIX}_extreme_F.txt"
echo "=================================================================="