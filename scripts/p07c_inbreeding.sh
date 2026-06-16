#!/bin/bash
#SBATCH --job-name=inbreeding
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/p07c_inbreeding_%j.log

# ============================================================
# Inbreeding coefficient estimation
# Run AFTER p07a (individual het array) and p07b (pop theta)
# are both complete.
#
# Fixed version:
#   - Base R only — no tidyverse dependency
#   - module load R for cluster compatibility
#   - NaN-safe IBS diagonal check (was crashing on NaN ibsMat)
#   - F_IBS skipped gracefully if matrix is unreliable
#   - Population summary written before method comparison
# ============================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
module load R/4.4.1-foss-2022b 2>/dev/null || true

set -euo pipefail

OUT_DIR=${OUT_DIR:-physalia/angsd_results}
METADATA="$OUT_DIR/metadata.tsv"
HET_DIR="$OUT_DIR/heterozygosity_corrected"
THETA_DIR="$OUT_DIR/theta_corrected"
INBREEDING_DIR="$OUT_DIR/inbreeding_corrected"
BAM_LIST="$OUT_DIR/all_samples.list"

mkdir -p "$INBREEDING_DIR"

echo "=================================================================="
echo "Inbreeding Coefficient Estimation"
echo "=================================================================="

Rscript --vanilla - \
    "$OUT_DIR" "$HET_DIR" "$THETA_DIR" "$METADATA" \
    "$INBREEDING_DIR" "$BAM_LIST" << 'REOF'

args          <- commandArgs(trailingOnly = TRUE)
out_dir       <- args[1]
het_dir       <- args[2]
theta_dir     <- args[3]
metadata_file <- args[4]
inbr_dir      <- args[5]
bam_list_file <- args[6]

metadata <- read.table(metadata_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)

# ============================================================
# 1. Load individual heterozygosity
# ============================================================
cat("\n=== Loading individual heterozygosity (1-sample SFS method) ===\n")

het_files <- list.files(het_dir, pattern = "\\.het$", full.names = TRUE)
cat("Found", length(het_files), "individual .het files\n")

if (length(het_files) == 0) {
    stop("ERROR: No .het files found in ", het_dir,
         "\nMake sure p07a array job completed successfully.")
}

het_list <- lapply(het_files, function(f) {
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
})
het_data <- do.call(rbind, het_list)
cat("Loaded heterozygosity for", nrow(het_data), "individuals\n")

het_data$sample_prefix <- sub("_.*", "", het_data$sample_id)
metadata$sample_prefix <- sub("_.*", "", metadata$sample_id)

het_data <- merge(het_data,
                  metadata[, c("sample_prefix", "population", "site", "year")],
                  by = "sample_prefix", all.x = TRUE)

het_data$pop_map <- paste0("pop_map_", het_data$site)

unmatched <- sum(is.na(het_data$site))
if (unmatched > 0) {
    cat("WARNING:", unmatched, "individuals could not be matched to metadata\n")
    print(het_data[is.na(het_data$site), "sample_id", drop = FALSE])
}

cat("\nIndividual heterozygosity by population:\n")
het_clean  <- het_data[!is.na(het_data$site), ]
pop_groups <- split(het_clean$heterozygosity, het_clean$pop_map)
het_pop_summary <- data.frame(
    pop_map  = names(pop_groups),
    n        = sapply(pop_groups, length),
    mean_het = round(sapply(pop_groups, mean,   na.rm = TRUE), 6),
    sd_het   = round(sapply(pop_groups, sd,     na.rm = TRUE), 6),
    min_het  = round(sapply(pop_groups, min,    na.rm = TRUE), 6),
    max_het  = round(sapply(pop_groups, max,    na.rm = TRUE), 6),
    stringsAsFactors = FALSE
)
print(het_pop_summary[order(het_pop_summary$pop_map), ])

write.table(het_data,
            file.path(inbr_dir, "individual_heterozygosity_corrected.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================
# 2. Load corrected population theta_pi
# ============================================================
cat("\n=== Loading corrected population theta_pi ===\n")

theta_file <- file.path(theta_dir, "theta_summary_corrected.tsv")
if (!file.exists(theta_file)) {
    stop("ERROR: Corrected theta summary not found: ", theta_file,
         "\nMake sure p07b completed successfully.")
}

theta_data <- read.table(theta_file, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
cat("Loaded theta for", nrow(theta_data), "populations\n")

# ============================================================
# METHOD 1: F_HET
# F = 1 - (H_observed / H_expected)
# H_obs  = individual heterozygosity from 1-sample SFS
# H_exp  = population theta_pi per site (corrected windowed thetaStat)
# ============================================================
cat("\n=== Method 1: F_HET (individual het vs population pi) ===\n")

theta_lookup <- setNames(theta_data$theta_pi_per_site,
                         theta_data$population)

f_het            <- het_clean
f_het$H_obs      <- f_het$heterozygosity
f_het$H_expected <- theta_lookup[f_het$pop_map]
f_het$F_HET      <- pmax(-1, pmin(1, 1 - (f_het$H_obs / f_het$H_expected)))

f_het <- f_het[order(-f_het$F_HET),
               c("sample_id", "pop_map", "site", "year",
                 "n_sites_total", "H_obs", "H_expected", "F_HET")]

cat("F_HET summary:\n")
cat("  Mean:  ", round(mean(f_het$F_HET,   na.rm = TRUE), 4), "\n")
cat("  Median:", round(median(f_het$F_HET, na.rm = TRUE), 4), "\n")
cat("  Range: ", round(min(f_het$F_HET,   na.rm = TRUE), 4),
    "to", round(max(f_het$F_HET, na.rm = TRUE), 4), "\n")
cat("  F > 0.10 (moderate inbreeding):",
    sum(f_het$F_HET > 0.10, na.rm = TRUE), "\n")
cat("  F > 0.25 (high inbreeding):",
    sum(f_het$F_HET > 0.25, na.rm = TRUE), "\n")

write.table(f_het,
            file.path(inbr_dir, "F_HET_individual_corrected.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================
# METHOD 2: F from IBS diagonal
# Only attempted if ibsMat exists and diagonal is valid (no NaN)
# The ibsMat for this dataset contains NaN — F_IBS is skipped
# ============================================================
cat("\n=== Method 2: F from IBS diagonal ===\n")

ibs_file <- file.path(out_dir, "pop_map_ALL.ibsMat")
f_ibs    <- NULL

if (!file.exists(ibs_file)) {
    cat("WARNING: IBS matrix not found — skipping F_IBS\n")
} else {
    ibs_mat   <- as.matrix(read.table(ibs_file, header = FALSE))
    diag_vals <- diag(ibs_mat)

    # NaN-safe check — must test for NaN before any comparison
    has_nan  <- any(is.nan(diag_vals))
    has_na   <- any(is.na(diag_vals))
    all_zero <- !has_nan && !has_na && all(diag_vals == 0)

    cat("IBS diagonal range:",
        round(min(diag_vals, na.rm = TRUE), 4), "to",
        round(max(diag_vals, na.rm = TRUE), 4), "\n")

    if (has_nan || has_na) {
        cat("WARNING: IBS diagonal contains NaN/NA values — matrix unreliable.",
            "Skipping F_IBS.\n")
        f_ibs <- NULL
    } else if (all_zero) {
        cat("WARNING: IBS diagonal is all zero — matrix invalid.",
            "Skipping F_IBS.\n")
        f_ibs <- NULL
    } else {
        if (mean(diag_vals) < 0.5) {
            cat("WARNING: IBS diagonal values suspiciously low (mean =",
                round(mean(diag_vals), 4), ") — interpret with caution\n")
        }

        bam_paths  <- readLines(bam_list_file)
        sample_ids <- sub("_merged.*", "", basename(bam_paths))
        n_use      <- min(length(sample_ids), length(diag_vals))

        f_ibs <- data.frame(
            sample_id    = sample_ids[1:n_use],
            IBS_diagonal = diag_vals[1:n_use],
            F_IBS        = diag_vals[1:n_use] - 1,
            stringsAsFactors = FALSE
        )

        f_ibs$sample_prefix <- sub("_.*", "", f_ibs$sample_id)
        f_ibs <- merge(f_ibs,
                       metadata[, c("sample_prefix", "population",
                                    "site", "year")],
                       by = "sample_prefix", all.x = TRUE)
        f_ibs <- f_ibs[order(-f_ibs$F_IBS), ]

        cat("F_IBS summary:\n")
        cat("  Mean: ", round(mean(f_ibs$F_IBS, na.rm = TRUE), 4), "\n")
        cat("  Range:", round(min(f_ibs$F_IBS,  na.rm = TRUE), 4),
            "to", round(max(f_ibs$F_IBS, na.rm = TRUE), 4), "\n")

        write.table(f_ibs,
                    file.path(inbr_dir, "F_IBS_individual_corrected.txt"),
                    sep = "\t", row.names = FALSE, quote = FALSE)
    }
}

# ============================================================
# Population-level summary
# Written regardless of whether F_IBS succeeded
# ============================================================
cat("\n=== Population-level inbreeding summary ===\n")

pop_split <- split(f_het, f_het$pop_map)

pop_summary <- do.call(rbind, lapply(names(pop_split), function(p) {
    d <- pop_split[[p]]
    data.frame(
        pop_map        = p,
        n_individuals  = nrow(d),
        mean_F_HET     = round(mean(d$F_HET,    na.rm = TRUE), 4),
        sd_F_HET       = round(sd(d$F_HET,      na.rm = TRUE), 4),
        median_F_HET   = round(median(d$F_HET,  na.rm = TRUE), 4),
        n_F_above_0.10 = sum(d$F_HET > 0.10,    na.rm = TRUE),
        n_F_above_0.25 = sum(d$F_HET > 0.25,    na.rm = TRUE),
        mean_H_obs     = round(mean(d$H_obs,     na.rm = TRUE), 6),
        mean_H_exp     = round(mean(d$H_expected,na.rm = TRUE), 6),
        stringsAsFactors = FALSE
    )
}))
pop_summary <- pop_summary[order(-pop_summary$mean_F_HET), ]

if (!is.null(f_ibs)) {
    ibs_split   <- split(f_ibs$F_IBS, paste0("pop_map_", f_ibs$site))
    ibs_summary <- data.frame(
        pop_map    = names(ibs_split),
        mean_F_IBS = round(sapply(ibs_split, mean, na.rm = TRUE), 4),
        sd_F_IBS   = round(sapply(ibs_split, sd,   na.rm = TRUE), 4),
        stringsAsFactors = FALSE
    )
    pop_summary <- merge(pop_summary, ibs_summary,
                         by = "pop_map", all.x = TRUE)
}

write.table(pop_summary,
            file.path(inbr_dir, "F_population_summary_corrected.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n")
print(pop_summary)

# ============================================================
# Method comparison (only if F_IBS succeeded)
# ============================================================
if (!is.null(f_ibs)) {
    cat("\n=== Method comparison ===\n")
    comp  <- merge(f_het[,  c("sample_id", "pop_map", "F_HET")],
                   f_ibs[, c("sample_id", "F_IBS")],
                   by = "sample_id")
    cor_r <- cor(comp$F_HET, comp$F_IBS, use = "complete.obs")
    cat("Pearson correlation F_HET vs F_IBS: r =", round(cor_r, 3), "\n")
    write.table(comp,
                file.path(inbr_dir, "F_method_comparison_corrected.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
}

# ============================================================
# Summary report
# ============================================================
report_lines <- c(
    "Inbreeding Coefficient Analysis — Corrected",
    "============================================",
    "",
    "Method 1: F_HET",
    "  H_observed  = individual heterozygosity from 1-sample SFS (realSFS -fold 1)",
    "  H_expected  = population theta_pi per site from corrected windowed thetaStat",
    "  F           = 1 - (H_obs / H_exp), clipped to [-1, 1]",
    "",
    "Method 2: F_IBS",
    "  Derived from IBS matrix diagonal (pop_map_ALL.ibsMat)",
    "  F = IBS_diagonal - 1",
    "  NOTE: Skipped for this dataset — ibsMat diagonal contains NaN values",
    "",
    sprintf("Individuals analysed: %d", nrow(f_het)),
    sprintf("F_HET mean:           %.4f", mean(f_het$F_HET, na.rm = TRUE)),
    sprintf("F_HET range:          %.4f to %.4f",
            min(f_het$F_HET, na.rm = TRUE),
            max(f_het$F_HET, na.rm = TRUE)),
    sprintf("F > 0.10 (moderate):  %d individuals",
            sum(f_het$F_HET > 0.10, na.rm = TRUE)),
    sprintf("F > 0.25 (high):      %d individuals",
            sum(f_het$F_HET > 0.25, na.rm = TRUE)),
    "",
    "Interpretation:",
    "  F < 0  : more heterozygous than population average (outbred)",
    "  F = 0  : average heterozygosity for population",
    "  F > 0.1: moderate inbreeding",
    "  F > 0.25: high inbreeding (equivalent to half-sibling mating)"
)

writeLines(report_lines,
           file.path(inbr_dir, "inbreeding_summary_corrected.txt"))

cat("\nAll results written to:", inbr_dir, "\n")
REOF

echo ""
echo "=================================================================="
echo "Inbreeding analysis complete"
echo "Results in: $INBREEDING_DIR"
echo "=================================================================="