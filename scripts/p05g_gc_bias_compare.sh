#!/bin/bash
#SBATCH --job-name=gc_bias_compare
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/gc_bias_compare_%j.log

# =============================================================================
# p05g_gc_bias_compare.sh — paired comparison of GC-content distribution
# (from samtools stats GCF/GCL lines) between amplified and non-amplified
# samples of the same individual.
#
# Computes, per sample: mean read GC%, and the proportion of reads in the
# extreme GC bins (<20% or >80%) — a standard way to summarise GC-bias
# dropout without needing a full distributional comparison.
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen   # R lives in this env
module load R/4.4.1-foss-2022b 2>/dev/null || true
set -eu

GC_DIR="physalia/angsd_results/amp_bias/gc_bias"

Rscript --vanilla - "$GC_DIR" << 'REOF'
args   <- commandArgs(trailingOnly = TRUE)
gc_dir <- args[1]

individuals <- c("NS40","SI7","NS14","NS41","SI18","CAI6","SI59",
                 "NS1","NS30","NS32")   # the 10 raw ("o") vs amplified pairs; T/aT are amplified too

summarise_gc <- function(label) {
    f <- file.path(gc_dir, paste0(label, ".samtools_stats.txt"))
    if (!file.exists(f)) { cat("WARNING: missing", f, "\n"); return(NULL) }
    lines <- readLines(f)
    gcf <- lines[grepl("^GCF", lines)]
    gcl <- lines[grepl("^GCL", lines)]
    parse_bins <- function(x) {
        parts <- do.call(rbind, strsplit(x, "\t"))
        data.frame(gc = as.numeric(parts[,2]), pct = as.numeric(parts[,3]))
    }
    df <- rbind(parse_bins(gcf), parse_bins(gcl))
    df <- aggregate(pct ~ gc, data = df, sum)
    df$pct <- df$pct / sum(df$pct) * 100   # renormalise combined F+L to 100%

    mean_gc    <- sum(df$gc * df$pct) / 100
    extreme_pct <- sum(df$pct[df$gc < 20 | df$gc > 80])
    data.frame(mean_gc = mean_gc, extreme_gc_pct = extreme_pct)
}

rows <- do.call(rbind, lapply(individuals, function(ind) {
    amp  <- summarise_gc(paste0(ind, "_amplified"))
    orig <- summarise_gc(paste0(ind, "_non_amplified"))
    if (is.null(amp) || is.null(orig)) return(NULL)
    data.frame(individual = ind,
               mean_gc_amplified  = amp$mean_gc,
               mean_gc_original   = orig$mean_gc,
               extreme_gc_pct_amplified = amp$extreme_gc_pct,
               extreme_gc_pct_original  = orig$extreme_gc_pct)
}))

cat("=== Paired GC-content summary (samtools stats) ===\n")
print(rows, digits = 4)

cat("\n=== Mean read GC%: paired Wilcoxon ===\n")
print(wilcox.test(rows$mean_gc_amplified, rows$mean_gc_original, paired = TRUE))
cat("Mean GC%: amplified =", round(mean(rows$mean_gc_amplified),3),
    " non_amplified =", round(mean(rows$mean_gc_original),3), "\n")

cat("\n=== Extreme-GC read proportion (<20% or >80% GC): paired Wilcoxon ===\n")
print(wilcox.test(rows$extreme_gc_pct_amplified, rows$extreme_gc_pct_original, paired = TRUE))
cat("Mean extreme-GC%: amplified =", round(mean(rows$extreme_gc_pct_amplified),3),
    " non_amplified =", round(mean(rows$extreme_gc_pct_original),3), "\n")

out_file <- file.path(gc_dir, "gc_bias_paired_summary.tsv")
write.table(rows, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nWritten to:", out_file, "\n")
REOF
