#!/bin/bash
#SBATCH --job-name=final_results_plots
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/final_results_plots_%j.log

# =============================================================================
# p_final_results_plots.sh — generates every plot NOT already produced by
# earlier scripts (Tajima's D, heterozygosity, PCAngsd F, PCA, admixture,
# sex ratio are already done): FST heatmap, ROH summary, HWE summary,
# lcMLkin relatedness summary, amplification-bias summary.
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
# R lives inside owenspopgen itself (installed post-maintenance, since the
# R/4.4.1-foss-2022b module no longer exists on the cluster) - no module load needed.
set -eu

OUT_DIR="physalia/angsd_results"
PLOT_DIR="$OUT_DIR/final_plots"
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

site_labels_stacked <- c(CAI = "Castle\nIsland", NS = "Nonsuch\nIsland",
                         SB = "Sinky\nBay",       SI = "Southampton\nIsland")

# ---------------- 1. FST heatmap ----------------
pairs <- list(c("CAI","NS"), c("CAI","SB"), c("CAI","SI"),
              c("NS","SB"), c("NS","SI"), c("SB","SI"))
fst_vals <- rbindlist(lapply(pairs, function(p) {
    f <- file.path(out_dir, "fst_corrected", paste0(p[1], "_", p[2], ".fst.global.txt"))
    if (!file.exists(f)) return(NULL)
    v <- as.numeric(strsplit(readLines(f), "\t")[[1]][2])
    data.table(pop1 = p[1], pop2 = p[2], fst = v)
}))
if (nrow(fst_vals) > 0) {
    # mirror for a symmetric matrix plot
    mirrored <- rbind(fst_vals, fst_vals[, .(pop1 = pop2, pop2 = pop1, fst)])
    pops <- c("CAI","NS","SB","SI")
    mirrored[, pop1 := factor(pop1, levels = pops)]
    mirrored[, pop2 := factor(pop2, levels = pops)]
    p <- ggplot(mirrored, aes(pop1, pop2, fill = fst)) +
        geom_tile(color = "white") +
        geom_text(aes(label = sprintf("%.3f", fst)), size = 5) +
        scale_fill_gradient(low = "white", high = "firebrick", name = expression(italic(F)[ST])) +
        scale_x_discrete(labels = site_labels_stacked) +
        scale_y_discrete(labels = site_labels_stacked) +
        labs(x = NULL, y = NULL) +
        theme_report(13) + coord_fixed()
    ggsave(file.path(plot_dir, "01_fst_heatmap.pdf"), p, width = 6.5, height = 6)
    cat("FST heatmap written\n")
} else cat("WARNING: no FST files found, skipping heatmap\n")

# ---------------- 2. ROH population summary ----------------
roh_f <- file.path(out_dir, "roh", "population", "population_ROH_summary.tsv")
if (file.exists(roh_f)) {
    roh <- fread(roh_f)
    p <- ggplot(roh, aes(x = population, y = mean_theta, fill = population)) +
        geom_col(alpha = 0.8) +
        geom_errorbar(aes(ymin = mean_theta - sd_theta, ymax = mean_theta + sd_theta), width = 0.2) +
        scale_x_discrete(labels = site_labels_stacked) +
        labs(x = NULL, y = expression(paste("Runs of Homozygosity (", theta, ")"))) +
        theme_report(12) + theme(legend.position = "none")
    ggsave(file.path(plot_dir, "02_roh_theta_by_population.pdf"), p, width = 7, height = 6)
    cat("ROH theta plot written\n")
} else cat("WARNING: ROH summary not found, skipping\n")

# ---------------- 3. HWE summary ----------------
hwe_f <- file.path(out_dir, "hwe", "hwe_summary_all_populations.tsv")
if (file.exists(hwe_f)) {
    hwe <- fread(hwe_f)
    p <- ggplot(hwe, aes(x = population, y = pct_sig_p005, fill = population)) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = 5, linetype = "dashed", colour = "grey30") +
        scale_x_discrete(labels = site_labels_stacked) +
        labs(x = NULL, y = "% of sites deviating significantly from HWE") +
        theme_report(12) + theme(legend.position = "none")
    ggsave(file.path(plot_dir, "03_hwe_deviation_by_population.pdf"), p, width = 7, height = 6)
    cat("HWE plot written\n")
} else cat("WARNING: HWE summary not found, skipping\n")

# ---------------- 4. lcMLkin relatedness summary ----------------
lck_f <- file.path(out_dir, "lcmlkin", "lcmlkin_all_populations.tsv")
if (file.exists(lck_f)) {
    lck <- fread(lck_f)
    classify <- function(Z0, r) {
        fifelse(r >= 0.45, "Duplicate",
        fifelse(Z0 <= 0.10 & r >= 0.40, "1st degree",
        fifelse(Z0 >= 0.15 & Z0 <= 0.35 & r >= 0.35, "Full sibling",
        fifelse(Z0 >= 0.40 & Z0 <= 0.60, "2nd degree",
        fifelse(Z0 >= 0.65 & Z0 <= 0.85, "3rd degree", "Unrelated")))))
    }
    lck[, relationship := classify(Z0ag, PI_HATag)]
    lck[, relationship := factor(relationship,
        levels = c("Duplicate","1st degree","Full sibling","2nd degree","3rd degree","Unrelated"))]

    summary_tab <- lck[relationship != "Unrelated", .N, by = .(population, relationship)]
    p <- ggplot(summary_tab, aes(x = population, y = N, fill = relationship)) +
        geom_col(position = "stack") +
        scale_x_discrete(labels = site_labels_stacked) +
        labs(x = NULL, y = "Number of pairs", fill = "Relationship") +
        theme_report(12)
    ggsave(file.path(plot_dir, "04_relatedness_by_population.pdf"), p, width = 7.5, height = 6)

    fwrite(lck[relationship != "Unrelated"][order(population, -PI_HATag)],
           file.path(plot_dir, "relatedness_close_pairs.tsv"), sep = "\t")
    cat("Relatedness plot written\n")
} else cat("WARNING: lcMLkin combined results not found, skipping\n")

# ---------------- 5. Amplification bias summary ----------------
amp_f <- file.path(out_dir, "amp_bias", "amplification_bias_data.tsv")
gc_f  <- file.path(out_dir, "amp_bias", "gc_bias", "gc_bias_paired_summary.tsv")
valid <- c("NS40","SI7","NS14","NS41","SI18","CAI6","SI59","NS1","NS30","NS32")
if (file.exists(amp_f) && file.exists(gc_f)) {
    amp <- fread(amp_f)[individual %in% valid]
    gc  <- fread(gc_f)[individual %in% valid]

    p1 <- ggplot(amp, aes(x = amplified, y = breadth, fill = amplified)) +
        geom_boxplot(alpha = 0.7) +
        geom_line(aes(group = individual), colour = "grey50", alpha = 0.5) +
        geom_point(size = 1.5) +
        labs(x = NULL, y = "genome coverage breadth") +
        theme_report(12) + theme(legend.position = "none")
    ggsave(file.path(plot_dir, "05_amp_bias_breadth.pdf"), p1, width = 6, height = 5.5)

    gc_long <- melt(gc[, .(individual, mean_gc_amplified, mean_gc_original)],
                    id.vars = "individual", variable.name = "group", value.name = "mean_gc")
    gc_long[, group := fifelse(group == "mean_gc_amplified", "amplified", "non_amplified")]
    p2 <- ggplot(gc_long, aes(x = group, y = mean_gc, fill = group)) +
        geom_boxplot(alpha = 0.7) +
        geom_line(aes(group = individual), colour = "grey50", alpha = 0.5) +
        geom_point(size = 1.5) +
        labs(x = NULL, y = "mean GC%") +
        theme_report(12) + theme(legend.position = "none")
    ggsave(file.path(plot_dir, "06_amp_bias_gc.pdf"), p2, width = 6, height = 5.5)
    cat("Amplification bias plots written\n")
} else cat("WARNING: amplification bias data not found, skipping\n")

cat("\nAll available plots written to:", plot_dir, "\n")
REOF

# -----------------------------------------------------------------------
# Export every PDF in $PLOT_DIR to PNG (300dpi) and SVG. pdftoppm/pdf2svg
# live in the owenspopgen env (already active here, but re-assert it).
# -----------------------------------------------------------------------
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
command -v pdftoppm >/dev/null 2>&1 || { echo "ERROR: pdftoppm not found in owenspopgen env"; exit 1; }
command -v pdf2svg  >/dev/null 2>&1 || { echo "ERROR: pdf2svg not found in owenspopgen env"; exit 1; }

for pdf in "$PLOT_DIR"/*.pdf; do
    [[ -f "$pdf" ]] || continue
    base="${pdf%.pdf}"
    pdftoppm -png -r 300 -singlefile "$pdf" "$base"
    pdf2svg "$pdf" "${base}.svg"
    echo "  exported: ${base}.png ${base}.svg"
done
