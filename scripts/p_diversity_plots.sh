#!/bin/bash
#SBATCH --job-name=diversity_plots
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/diversity_plots_%j.log

# =============================================================================
# p_diversity_plots.sh — population-level summary plots for heterozygosity
# and inbreeding, for the team report. Site plots show the 2023 and 2024
# sequencing runs side by side; plot 05 is 2024-only theta by site. Reads only already-computed outputs
# (p07a .het files, p07c F_HET table, PCAngsd individual F) — no new ANGSD.
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
# R lives inside owenspopgen itself (installed post-maintenance, since the
# R/4.4.1-foss-2022b module no longer exists on the cluster) - no module load needed.
set -eu

OUT_DIR=${OUT_DIR:-physalia/angsd_results}
PLOT_DIR="$OUT_DIR/diversity_plots"
mkdir -p "$PLOT_DIR"

Rscript --vanilla - "$OUT_DIR" "$PLOT_DIR" << 'REOF'
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

theme_report <- function(base_size = 13) {
  theme_classic(base_size = base_size) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key = element_rect(fill = "white", colour = NA))
}

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]; plot_dir <- args[2]

# stacked (two-line) full site names for x-axis tick labels
site_labels_stacked <- c(CAI = "Castle\nIsland", NS = "Nonsuch\nIsland",
                         SB = "Sinky\nBay",       SI = "Southampton\nIsland")

meta <- fread(file.path(out_dir, "metadata.tsv"))
meta[, sample_prefix := sub("_.*", "", sample_id)]
main_pops <- c("CAI", "NS", "SB", "SI")   # match theta/report convention; small pops shown but flagged

# Sequencing run: samples numbered 1-63 were sequenced in the 2024 run, all others (and every
# library ending in "o") in the 2023 run. Run shifts individual heterozygosity and F, so every
# site plot shows the two runs side by side instead of pooling them.
seq_run <- function(prefix) {
  num <- suppressWarnings(as.integer(sub("[^0-9].*", "", prefix)))
  factor(ifelse(!grepl("o$", prefix) & !is.na(num) & num <= 63, "2024", "2023"), levels = c("2023", "2024"))
}
run_cols <- c("2023" = "#9E9E9E", "2024" = "#2C7FB8")
run_box <- function(dt, yvar, ylab, hline = NULL) {
  p <- ggplot(dt, aes(x = site, y = .data[[yvar]], fill = run)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8, preserve = "single"), width = 0.7) +
    geom_point(aes(colour = run), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               alpha = 0.5, size = 1) +
    scale_fill_manual(values = run_cols, name = "Sequencing run", drop = FALSE) +
    scale_colour_manual(values = run_cols, guide = "none", drop = FALSE) +
    scale_x_discrete(labels = site_labels_stacked) +
    labs(x = NULL, y = ylab) +
    theme_report(12) + theme(legend.position = "top")
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", colour = "grey40")
  p
}

# ---------------- 1. Genome-wide individual heterozygosity ----------------
het_files <- list.files(file.path(out_dir, "heterozygosity_corrected"), pattern = "\\.het$", full.names = TRUE)
het <- rbindlist(lapply(het_files, function(f) tryCatch(fread(f), error=function(e) NULL)))
setnames(het, "sample_id", "sample", skip_absent = TRUE)
het[, sample_prefix := sub("_.*", "", sample)]
het <- merge(het, meta[, .(sample_prefix, site)], by = "sample_prefix", all.x = TRUE)
het[, is_main := site %in% main_pops]
het[, run := seq_run(sample_prefix)]

p1 <- run_box(het[is_main == TRUE], "heterozygosity", "Heterozygosity (het sites / total sites)")
ggsave(file.path(plot_dir, "01_heterozygosity_by_population.pdf"), p1, width = 7, height = 5.5)

p1b <- ggplot(het[is_main == TRUE], aes(x = heterozygosity, fill = run)) +
    geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
    scale_fill_manual(values = run_cols, name = "Sequencing run", drop = FALSE) +
    facet_wrap(~site, scales = "free_y", ncol = 2,
               labeller = labeller(site = site_labels_stacked)) +
    labs(x = "Heterozygosity", y = "count") +
    theme_report(11) + theme(legend.position = "top")
ggsave(file.path(plot_dir, "02_heterozygosity_distribution.pdf"), p1b, width = 8, height = 6)

cat("=== Heterozygosity summary by population ===\n")
print(het[is_main == TRUE, .(n = .N, mean_het = mean(heterozygosity, na.rm=TRUE),
                              sd_het = sd(heterozygosity, na.rm=TRUE)), by = .(site, run)][order(site, run)])

# ---------------- 2. Individual F: PCAngsd (preferred; structure-corrected) ----------------
pcangsd_f_file <- file.path(out_dir, "pcangsd_inbreed_ind.inbreed.samples")
pcangsd_list   <- file.path(out_dir, "pop_map_lists", "pop_map_ALL.list")
if (file.exists(pcangsd_f_file) && file.exists(pcangsd_list)) {
    samples <- sub("_merged.*", "", basename(readLines(pcangsd_list)))
    pf <- fread(pcangsd_f_file, header = FALSE)
    if (nrow(pf) == length(samples)) {
        pf[, sample := samples]
        pf[, sample_prefix := sub("_.*", "", sample)]
        pf <- merge(pf, meta[, .(sample_prefix, site)], by = "sample_prefix", all.x = TRUE)
        pf[, is_main := site %in% main_pops]
        pf[, run := seq_run(sample_prefix)]

        p2 <- run_box(pf[is_main == TRUE], "V1",
                      expression(paste("PCAngsd inbreeding coefficient (", italic(F), ")")), hline = 0)
        ggsave(file.path(plot_dir, "03_pcangsd_F_by_population.pdf"), p2, width = 7, height = 5.5)

        cat("\n=== PCAngsd individual F summary by population (PRIMARY inbreeding measure) ===\n")
        print(pf[is_main == TRUE, .(n = .N, mean_F = mean(V1, na.rm=TRUE),
                                     median_F = median(V1, na.rm=TRUE)), by = .(site, run)][order(site, run)])
        fwrite(pf, file.path(plot_dir, "pcangsd_F_by_individual.tsv"), sep = "\t")
    } else {
        cat("WARNING: PCAngsd F sample count mismatch — skipping\n")
    }
} else {
    cat("WARNING: PCAngsd inbreed-samples output not found — skipping (see notes: run pcangsd --inbreed-samples first)\n")
}

# ---------------- 3. F_HET (secondary; NOTE known upward-bias caveat) ----------------
fhet_file <- file.path(out_dir, "inbreeding_corrected", "F_HET_individual_corrected.txt")
if (file.exists(fhet_file)) {
    fh <- fread(fhet_file)
    fh <- fh[!is.na(F_HET) & site %in% main_pops]
    pref <- if ("sample_prefix" %in% names(fh)) fh$sample_prefix else sub("_.*", "", fh[[1]])
    fh[, run := seq_run(pref)]
    p3 <- run_box(fh, "F_HET", expression(italic(F)[HET]))
    ggsave(file.path(plot_dir, "04_F_HET_by_population_CAVEAT.pdf"), p3, width = 7, height = 5.5)
    cat("\nF_HET plot written (flagged with bias caveat in subtitle)\n")
}

# ---------------- 4. Population diversity by site, 2024 run only ----------------
# From bx_theta_by_year.sh (SITES="CAI NS SB SI" YEARS="2024" N=12 ... --collate).
th_file <- file.path(out_dir, "theta_by_year", "theta_by_site_year_2024_n12.tsv")
if (file.exists(th_file)) {
    th <- fread(th_file)[site %in% main_pops]
    tl <- rbind(th[, .(site, stat = "theta[W]", value = thetaW_per_site)],
                th[, .(site, stat = "pi", value = pi_per_site)])
    tl[, stat := factor(stat, levels = c("theta[W]", "pi"))]
    p4 <- ggplot(tl, aes(x = site, y = value * 1e4, fill = stat)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
        geom_text(data = th, aes(x = site, y = pmax(thetaW_per_site, pi_per_site) * 1e4,
                                 label = sprintf("D = %.2f", mean_tajimaD)),
                  inherit.aes = FALSE, vjust = -0.6, size = 3.5) +
        scale_fill_manual(values = c("theta[W]" = "#2C7FB8", "pi" = "#7FCDBB"),
                          labels = c("theta[W]" = expression(theta[W]), "pi" = expression(pi)), name = NULL) +
        scale_x_discrete(labels = site_labels_stacked) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
        labs(x = NULL, y = expression(paste("Diversity per site (", "" %*% 10^-4, ")"))) +
        theme_report(12) + theme(legend.position = "top")
    ggsave(file.path(plot_dir, "05_theta_by_site_2024.pdf"), p4, width = 7, height = 5.5)
    cat("\n2024-only theta plot written\n")
} else {
    cat("\nNOTE:", th_file, "not found - skipping the 2024-only theta plot\n")
}

cat("\nAll diversity plots written to:", plot_dir, "\n")
REOF

# -----------------------------------------------------------------------
# Export every PDF in $PLOT_DIR to PNG (300dpi) and SVG.
# pdftoppm/pdf2svg live in the owenspopgen env - re-activate explicitly
# here since a caller may have activated a different env for the R step.
# -----------------------------------------------------------------------
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
command -v pdftoppm >/dev/null 2>&1 || { echo "ERROR: pdftoppm not found in owenspopgen env (conda install -c conda-forge poppler)"; exit 1; }
command -v pdf2svg  >/dev/null 2>&1 || { echo "ERROR: pdf2svg not found in owenspopgen env (conda install -c conda-forge pdf2svg)"; exit 1; }

for pdf in "$PLOT_DIR"/*.pdf; do
    [[ -f "$pdf" ]] || continue
    base="${pdf%.pdf}"
    pdftoppm -png -r 300 -singlefile "$pdf" "$base"
    pdf2svg "$pdf" "${base}.svg"
    echo "  exported: ${base}.png ${base}.svg"
done

echo ""
echo "Done. Plots (pdf/png/svg) in: $PLOT_DIR"
