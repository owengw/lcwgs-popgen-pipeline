#!/bin/bash
#SBATCH --job-name=inbreeding
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/p07c_inbreeding_%j.log

# =============================================================================
# p07c_inbreeding.sh — RagTag rewrite. F_HET = 1 - (H_obs/H_exp) per individual.
# H_obs from p07a's .het files (all-sites, 180 individuals).
# H_exp from p06's theta_pi (via the p06_theta_to_corrected.sh adapter).
# No EXCLUDE regex — reads sample identity from metadata.tsv, which was
# itself generated from the correct 180-sample bam list this session.
# Base R only, matches repo convention.
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
module load R/4.4.1-foss-2022b 2>/dev/null || true
set -eu

OUT_DIR=${OUT_DIR:-physalia/angsd_results}
METADATA="$OUT_DIR/metadata.tsv"
HET_DIR="$OUT_DIR/heterozygosity_corrected"
THETA_DIR="$OUT_DIR/theta_corrected"
INBREEDING_DIR="$OUT_DIR/inbreeding_corrected"
IBS_FILE="$OUT_DIR/pop_map_ALL.ibsMat"   # optional; produced by p06c2_pca.sh if run

mkdir -p "$INBREEDING_DIR"

echo "=================================================================="
echo "Inbreeding Coefficient Estimation (F_HET, all-sites)"
echo "=================================================================="

Rscript --vanilla - "$OUT_DIR" "$HET_DIR" "$THETA_DIR" "$METADATA" "$INBREEDING_DIR" "$IBS_FILE" << 'REOF'
args          <- commandArgs(trailingOnly = TRUE)
out_dir       <- args[1]
het_dir       <- args[2]
theta_dir     <- args[3]
metadata_file <- args[4]
inbr_dir      <- args[5]
ibs_file      <- args[6]

metadata <- read.table(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# ============================================================
# 1. Individual heterozygosity (from p07a)
# ============================================================
cat("\n=== Loading individual heterozygosity ===\n")
het_files <- list.files(het_dir, pattern = "\\.het$", full.names = TRUE)
cat("Found", length(het_files), "individual .het files\n")
if (length(het_files) == 0) stop("No .het files in ", het_dir, " — run p07a first")

het_data <- do.call(rbind, lapply(het_files, function(f)
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)))
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

# ============================================================
# 2. Population theta_pi (from p06, via adapter)
# ============================================================
cat("\n=== Loading population theta_pi ===\n")
theta_file <- file.path(theta_dir, "theta_summary_corrected.tsv")
if (!file.exists(theta_file)) {
    stop("theta_summary_corrected.tsv not found: ", theta_file,
         "\nRun p06_theta_to_corrected.sh first")
}
theta_data <- read.table(theta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Loaded theta for", nrow(theta_data), "populations\n")

# ============================================================
# METHOD 1: F_HET
# ============================================================
cat("\n=== F_HET (individual het vs population theta_pi) ===\n")

theta_lookup <- setNames(theta_data$theta_pi_per_site, theta_data$population)

f_het <- het_data[!is.na(het_data$site), ]
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
cat("  F > 0.10 (moderate inbreeding):", sum(f_het$F_HET > 0.10, na.rm = TRUE), "\n")
cat("  F > 0.25 (high inbreeding):",     sum(f_het$F_HET > 0.25, na.rm = TRUE), "\n")

write.table(f_het, file.path(inbr_dir, "F_HET_individual_corrected.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ============================================================
# METHOD 2: F_IBS (optional — only if the IBS matrix exists)
# ============================================================
cat("\n=== F_IBS (IBS matrix diagonal, if available) ===\n")
f_ibs <- NULL
if (nchar(ibs_file) > 0 && file.exists(ibs_file)) {
    ibs_mat   <- as.matrix(read.table(ibs_file, header = FALSE))
    diag_vals <- diag(ibs_mat)
    has_nan <- any(is.nan(diag_vals)); has_na <- any(is.na(diag_vals))
    if (has_nan || has_na) {
        cat("WARNING: IBS diagonal has NaN/NA — skipping F_IBS\n")
    } else {
        bam_list <- file.path(out_dir, "pop_map_lists", "pop_map_ALL.list")
        if (file.exists(bam_list)) {
            bam_paths  <- readLines(bam_list)
            sample_ids <- sub("_merged.*", "", basename(bam_paths))
            n_use <- min(length(sample_ids), length(diag_vals))
            f_ibs <- data.frame(sample_id = sample_ids[1:n_use],
                                 IBS_diagonal = diag_vals[1:n_use],
                                 F_IBS = diag_vals[1:n_use] - 1,
                                 stringsAsFactors = FALSE)
            f_ibs$sample_prefix <- sub("_.*", "", f_ibs$sample_id)
            f_ibs <- merge(f_ibs, metadata[, c("sample_prefix","population","site","year")],
                           by = "sample_prefix", all.x = TRUE)
            f_ibs <- f_ibs[order(-f_ibs$F_IBS), ]
            write.table(f_ibs, file.path(inbr_dir, "F_IBS_individual_corrected.txt"),
                        sep = "\t", row.names = FALSE, quote = FALSE)
            cat("F_IBS written (", nrow(f_ibs), "individuals )\n")
        } else {
            cat("pop_map_ALL.list not found — skipping F_IBS sample matching\n")
        }
    }
} else {
    cat("No IBS matrix at", ibs_file, "— skipping F_IBS (run p06c2_pca.sh if wanted)\n")
}

# ============================================================
# Population-level summary
# ============================================================
cat("\n=== Population-level inbreeding summary ===\n")
pop_split <- split(f_het, f_het$pop_map)
pop_summary <- do.call(rbind, lapply(names(pop_split), function(p) {
    d <- pop_split[[p]]
    data.frame(pop_map = p, n_individuals = nrow(d),
               mean_F_HET = round(mean(d$F_HET, na.rm=TRUE), 4),
               sd_F_HET   = round(sd(d$F_HET,   na.rm=TRUE), 4),
               median_F_HET = round(median(d$F_HET, na.rm=TRUE), 4),
               n_F_above_0.10 = sum(d$F_HET > 0.10, na.rm=TRUE),
               n_F_above_0.25 = sum(d$F_HET > 0.25, na.rm=TRUE),
               stringsAsFactors = FALSE)
}))
pop_summary <- pop_summary[order(-pop_summary$mean_F_HET), ]
write.table(pop_summary, file.path(inbr_dir, "F_population_summary_corrected.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
print(pop_summary)

cat("\nAll results written to:", inbr_dir, "\n")
REOF

echo ""
echo "=================================================================="
echo "Inbreeding analysis complete — results in $OUT_DIR/inbreeding_corrected"
echo "=================================================================="
