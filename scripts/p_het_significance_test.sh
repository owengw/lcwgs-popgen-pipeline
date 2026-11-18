#!/bin/bash
#SBATCH --job-name=het_sig_test
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/het_sig_test_%j.log

# =============================================================================
# p_het_significance_test.sh — tests whether individual genome-wide
# heterozygosity differs significantly among sampling sites.
# Kruskal-Wallis (all sites) + pairwise Wilcoxon post-hoc (BH-corrected),
# using the corrected p07a heterozygosity_corrected/*.het files.
# =============================================================================

source ~/.bash_profile
module load R/4.4.1-foss-2022b 2>/dev/null || true
export R_LIBS_USER=/mnt/parscratch/users/bi4og/R_libs
set -eu

OUT_DIR="physalia/angsd_results"

Rscript --vanilla - "$OUT_DIR" << 'REOF'
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]

het_files <- list.files(file.path(out_dir, "heterozygosity_corrected"), pattern = "\\.het$", full.names = TRUE)
het <- rbindlist(lapply(het_files, function(f) tryCatch(fread(f), error = function(e) NULL)))
setnames(het, "sample_id", "sample", skip_absent = TRUE)

meta <- fread(file.path(out_dir, "metadata.tsv"))
het[, sample_prefix := sub("_.*", "", sample)]
meta[, sample_prefix := sub("_.*", "", sample_id)]
het <- merge(het, meta[, .(sample_prefix, site)], by = "sample_prefix", all.x = TRUE)
het <- het[!is.na(heterozygosity) & !is.na(site)]

# only sites with n>=3 (matches convention used elsewhere in this pipeline)
n_per_site <- het[, .N, by = site]
keep_sites <- n_per_site[N >= 3, site]
het_f <- het[site %in% keep_sites]

cat("=== Sample sizes per site (n>=3 retained) ===\n")
print(n_per_site[order(-N)])

cat("\n=== Summary statistics per site ===\n")
summary_tab <- het_f[, .(n = .N, median_het = median(heterozygosity),
                          mean_het = mean(heterozygosity), sd_het = sd(heterozygosity),
                          se_het = sd(heterozygosity)/sqrt(.N)), by = site][order(-median_het)]
print(summary_tab)

cat("\n=== Kruskal-Wallis test: heterozygosity ~ site ===\n")
kw <- kruskal.test(heterozygosity ~ site, data = het_f)
print(kw)

if (kw$p.value < 0.05) {
    cat("\n=== Significant KW result — pairwise Wilcoxon post-hoc (BH-corrected) ===\n")
    ph <- pairwise.wilcox.test(het_f$heterozygosity, het_f$site, p.adjust.method = "BH")
    print(ph)
} else {
    cat("\nKW not significant at p<0.05 — no post-hoc comparisons performed.\n")
}

out_file <- file.path(out_dir, "diversity_plots", "heterozygosity_by_site_summary.tsv")
fwrite(summary_tab, out_file, sep = "\t")
cat("\nSummary written to:", out_file, "\n")
REOF
