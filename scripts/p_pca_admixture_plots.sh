#!/bin/bash
#SBATCH --job-name=pca_admix_plots
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/pca_admix_plots_%j.log

# =============================================================================
# p_pca_admixture_plots.sh — eigen-decomposes the PCAngsd covariance matrix
# and plots PC1/PC2/PC3 coloured by population, plots NGSadmix log-likelihood
# vs K (for choosing K), and admixture barplots for a few K values.
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
module load R/4.4.1-foss-2022b 2>/dev/null || true
export R_LIBS_USER=/mnt/parscratch/users/bi4og/R_libs
set -eu

OUT_DIR=${OUT_DIR:-physalia/angsd_results}
PLOT_DIR="$OUT_DIR/pca_admixture_plots"
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

# ---------------- PCA ----------------
cov <- as.matrix(read.table(file.path(out_dir, "pcangsd_ALL.cov")))
samples <- sub("_merged.*", "", basename(readLines(file.path(out_dir, "pop_map_lists/pop_map_ALL.list"))))
rownames(cov) <- colnames(cov) <- samples

meta <- fread(file.path(out_dir, "metadata.tsv"))
meta[, sample_prefix := sub("_.*", "", sample_id)]

eig <- eigen(cov)
var_pct <- eig$values / sum(eig$values) * 100
pcs <- data.table(sample = samples,
                  PC1 = eig$vectors[,1], PC2 = eig$vectors[,2],
                  PC3 = eig$vectors[,3], PC4 = eig$vectors[,4])
pcs[, sample_prefix := sub("_.*", "", sample)]
pcs <- merge(pcs, meta[, .(sample_prefix, population, site)], by = "sample_prefix", all.x = TRUE)

cat("Variance explained: PC1=", round(var_pct[1],1), "%  PC2=", round(var_pct[2],1),
    "%  PC3=", round(var_pct[3],1), "%\n")

p1 <- ggplot(pcs, aes(PC1, PC2, colour = site)) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2]),
         colour = "Site") +
    theme_report(12) + theme(legend.position = "right")
ggsave(file.path(plot_dir, "PCA_PC1_PC2.pdf"), p1, width = 8, height = 6)

p2 <- ggplot(pcs, aes(PC1, PC3, colour = site)) +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC3 (%.1f%%)", var_pct[3]),
         colour = "Site") +
    theme_report(12) + theme(legend.position = "right")
ggsave(file.path(plot_dir, "PCA_PC1_PC3.pdf"), p2, width = 8, height = 6)

fwrite(pcs, file.path(plot_dir, "PCA_coordinates.tsv"), sep = "\t")
cat("PCA plots + coordinates written\n")

# ---------------- PCAngsd admixture (--admix, preferred for lcWGS) ----------------
# NOTE: this reads pcangsd's own --admix output (K*.admix.Q), NOT NGSadmix's
# .qopt — PCAngsd's GL-based NMF estimator avoids NGSadmix's reliance on
# hard-called-genotype-adjacent likelihoods, which is more appropriate at
# lcWGS coverage. Run the pcangsd --admix jobs first (see launch notes).
admix_dir <- file.path(out_dir, "pcangsd_admix")
plot_pcangsd_admix_K <- function(K) {
    qfile <- file.path(admix_dir, sprintf("K%d.admix.Q", K))
    if (!file.exists(qfile)) { cat("No admix.Q for K=", K, "(expected", qfile, ")\n"); return(invisible(NULL)) }

    q <- fread(qfile, header = FALSE)
    q[, sample := samples[1:nrow(q)]]
    q[, sample_prefix := sub("_.*", "", sample)]
    q <- merge(q, meta[, .(sample_prefix, site)], by = "sample_prefix", all.x = TRUE)
    q <- q[order(site)]
    q[, sample := factor(sample, levels = sample)]

    long <- melt(q, id.vars = c("sample", "site"), measure.vars = paste0("V", 1:K),
                variable.name = "cluster", value.name = "prop")

    p <- ggplot(long, aes(x = sample, y = prop, fill = cluster)) +
        geom_col(width = 1) +
        facet_grid(~site, scales = "free_x", space = "free_x") +
        labs(x = NULL, y = "ancestry proportion") +
        theme_report(10) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              panel.spacing = unit(0.1, "lines"), legend.position = "none")
    ggsave(file.path(plot_dir, sprintf("PCAngsd_admix_K%d_barplot.pdf", K)), p, width = 12, height = 4)
    cat("PCAngsd admixture barplot K=", K, "written\n")
}

for (k in c(2, 3, 4)) plot_pcangsd_admix_K(k)

cat("\nAll PCA/admixture outputs written to:", plot_dir, "\n")
REOF

echo ""
echo "Done. Plots + tables in: $PLOT_DIR"
