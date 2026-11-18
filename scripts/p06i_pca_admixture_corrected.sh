#!/bin/bash
#SBATCH --job-name=pca_admixture_corrected
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/pca_admixture_corrected_%j.log

# =============================================================================
# p06i_pca_admixture_corrected.sh — PCA and admixture, rerun on the
# BAQ/C-corrected LD-pruned beagle (pcangsd_ALL_corrected.beagle.gz).
#
# PCA and admixture proportions are two SEPARATE PCAngsd runs:
#   - Default mode (no --admix) -> covariance matrix -> eigendecomposition
#     in R -> PC1/PC2/etc for PCA
#   - --admix mode, K=2-4 -> NMF-based admixture proportions
# Individual inbreeding (--inbreed-samples) also rerun here, since it
# depends on the same corrected beagle.
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/pcangsd310
set -eu

OUT_DIR="${OUT_DIR:-physalia/angsd_results}"
BEAGLE="$OUT_DIR/pcangsd_ALL_corrected.beagle.gz"
PLOT_DIR="$OUT_DIR/pca_admixture_plots_corrected"
PCANGSD=/mnt/parscratch/users/bi4og/conda_envs/pcangsd310/bin/pcangsd
mkdir -p "$PLOT_DIR"

[[ -f "$BEAGLE" ]] || { echo "ERROR: $BEAGLE not found"; exit 1; }

echo "=================================================================="
echo "1. PCA (default PCAngsd mode -> covariance matrix)"
echo "=================================================================="
"$PCANGSD" --beagle "$BEAGLE" --threads 8 \
    --out "$OUT_DIR/pcangsd_ALL_corrected"

echo "=================================================================="
echo "2. Admixture proportions, K=2-4"
echo "=================================================================="
mkdir -p "$OUT_DIR/pcangsd_admix_corrected"
for K in 2 3 4; do
    "$PCANGSD" --beagle "$BEAGLE" --threads 8 --admix --admix-K "$K" --admix-seed 42 \
        --out "$OUT_DIR/pcangsd_admix_corrected/K${K}"
done

echo "=================================================================="
echo "3. Individual inbreeding coefficients (structure-corrected)"
echo "=================================================================="
"$PCANGSD" --beagle "$BEAGLE" --threads 8 --inbreed-samples \
    --out "$OUT_DIR/pcangsd_inbreed_ind_corrected"

echo "=================================================================="
echo "4. PCA eigendecomposition, admixture summary, and plots"
echo "=================================================================="
# pcangsd (used above) has finished, so switch to owenspopgen for R - it has
# R installed directly (post-maintenance the R/4.4.1-foss-2022b module no
# longer exists, and pcangsd310 never had R anyway).
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

Rscript --vanilla - "$OUT_DIR" "$PLOT_DIR" << 'REOF'
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggh4x) })

theme_report <- function(base_size = 13) {
  theme_classic(base_size = base_size) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key = element_rect(fill = "white", colour = NA))
}

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]; plot_dir <- args[2]

# full site names (one line - used for legend keys, not axis ticks)
site_labels_full <- c(CAI = "Castle Island", NS = "Nonsuch Island",
                      SB = "Sinky Bay",       SI = "Southampton Island",
                      COI = "Cooper's Island", DH = "Daniel's Head",
                      SP = "Spittal Pond")

# ---- sample order (matches beagle column order = pop_map_ALL.list order) ----
samples <- sub("_merged.*", "", basename(readLines(file.path(out_dir, "pop_map_lists/pop_map_ALL.list"))))
meta <- fread(file.path(out_dir, "metadata.tsv"))
meta[, sample_prefix := sub("_.*", "", sample_id)]
samp_df <- data.table(sample = samples, sample_prefix = sub("_.*", "", samples))
samp_df <- merge(samp_df, meta[, .(sample_prefix, site)], by = "sample_prefix", all.x = TRUE, sort = FALSE)

# ---- 1. PCA from covariance matrix ----
cov <- as.matrix(fread(file.path(out_dir, "pcangsd_ALL_corrected.cov")))
e <- eigen(cov)
var_pct <- round(100 * e$values / sum(e$values), 1)

pcs <- data.table(sample_prefix = samp_df$sample_prefix, site = samp_df$site,
                  PC1 = e$vectors[,1], PC2 = e$vectors[,2], PC3 = e$vectors[,3])
fwrite(pcs, file.path(out_dir, "PCA_coordinates_corrected.tsv"), sep = "\t")
cat("PC1 variance explained:", var_pct[1], "%\n")
cat("PC2 variance explained:", var_pct[2], "%\n")
cat("PC3 variance explained:", var_pct[3], "%\n")

p_pc12 <- ggplot(pcs, aes(PC1, PC2, colour = site)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_colour_discrete(labels = site_labels_full) +
    labs(x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2]), colour = "Site") +
    theme_report(13)
ggsave(file.path(plot_dir, "PCA_PC1_PC2_corrected.pdf"), p_pc12, width = 7.5, height = 5.5)

p_pc13 <- ggplot(pcs, aes(PC1, PC3, colour = site)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_colour_discrete(labels = site_labels_full) +
    labs(x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC3 (%.1f%%)", var_pct[3]), colour = "Site") +
    theme_report(13)
ggsave(file.path(plot_dir, "PCA_PC1_PC3_corrected.pdf"), p_pc13, width = 7.5, height = 5.5)

# ---- 2. Admixture barplots + cluster summary, K=2-4 ----
for (K in 2:4) {
    qf <- file.path(out_dir, "pcangsd_admix_corrected", sprintf("K%d.admix.%d.Q", K, K))
    q <- fread(qf, header = FALSE)
    q[, sample_prefix := samp_df$sample_prefix]
    q <- merge(q, samp_df[, .(sample_prefix, site)], by = "sample_prefix", sort = FALSE)

    qcols <- paste0("V", 1:K)
    q[, dominant := qcols[apply(.SD, 1, which.max)], .SDcols = qcols]
    q[, max_prop := apply(.SD, 1, max), .SDcols = qcols]
    q[, admixed := max_prop < 0.9]

    cat("\n=== K =", K, "===\n")
    cat("Dominant cluster by site:\n")
    print(table(q$site, q$dominant))
    cat("Admixed individuals (max prop < 0.9):", sum(q$admixed), "/", nrow(q), "\n")

    q_long <- melt(q, id.vars = c("sample_prefix", "site"), measure.vars = qcols,
                  variable.name = "cluster", value.name = "proportion")
    q_long[, sample_prefix := factor(sample_prefix, levels = q[order(site, -get(qcols[1]))]$sample_prefix)]
    q_long[, cluster := factor(cluster, labels = paste("Cluster", seq_len(K)))]

    # facet strips carry the site name AND its sample size. Small-n sites
    # (Cooper's Island / Daniel's Head / Spittal Pond) have narrow facets, so
    # they get abbreviations instead of the full name
    site_labels_strip <- site_labels_full
    site_labels_strip[c("COI", "DH", "SP")] <- c("COI", "DH", "SP")
    n_by_site <- q[, .N, by = site]
    strip_labels <- setNames(
        sprintf("%s\n(n=%d)", site_labels_strip[n_by_site$site], n_by_site$N),
        n_by_site$site)
    # minimum panel width so small-n sites have room for their strip label;
    # site_order matches the alphabetical order facet_grid uses
    site_order <- levels(factor(q_long$site))
    panel_w <- pmax(n_by_site$N[match(site_order, n_by_site$site)], 6)

    p_admix <- ggplot(q_long, aes(x = sample_prefix, y = proportion, fill = cluster)) +
        geom_col(width = 1) +
        facet_grid(~site, scales = "free_x", space = "free_x",
                   labeller = labeller(site = strip_labels)) +
        labs(x = NULL, y = "Ancestry proportion") +
        theme_report(11) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              strip.text = element_text(size = 7.5),
              legend.position = "none") +
              force_panelsizes(cols = panel_w)
    ggsave(file.path(plot_dir, sprintf("PCAngsd_admix_K%d_barplot_corrected.pdf", K)), p_admix, width = 12, height = 4.5)

    fwrite(q, file.path(out_dir, sprintf("admixture_K%d_corrected_summary.tsv", K)), sep = "\t")
}

cat("\nAll corrected PCA/admixture outputs written to:", plot_dir, "\n")
REOF

# -----------------------------------------------------------------------
# Export every PDF in $PLOT_DIR to PNG (300dpi) and SVG.
# pdftoppm/pdf2svg live in owenspopgen, not pcangsd310 (active above) -
# switch envs explicitly for this step.
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
