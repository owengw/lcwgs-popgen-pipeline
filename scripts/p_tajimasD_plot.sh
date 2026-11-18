#!/bin/bash
#SBATCH --job-name=tajimasD_plot
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/tajimasD_plot_%j.log

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
# R lives inside owenspopgen itself (installed post-maintenance, since the
# R/4.4.1-foss-2022b module no longer exists on the cluster) - no module load needed.
set -eu

OUT_DIR=${OUT_DIR:-physalia/angsd_results}
PLOT_DIR="$OUT_DIR/diversity_plots"
mkdir -p "$PLOT_DIR"

Rscript --vanilla - "$OUT_DIR" "$PLOT_DIR" << 'REOF'
suppressPackageStartupMessages(library(ggplot2))

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

# use the sample-size-corrected (downsampled n=16) values — PRIMARY, per the
# earlier finding that raw full-n Tajima's D was inflated by sample size
f <- file.path(out_dir, "theta_downsampled", "theta_downsampled_summary.tsv")
dat <- read.delim(f, stringsAsFactors = FALSE)
dat <- dat[dat$analysis == "downsampled_n16", ]

p <- ggplot(dat, aes(x = population, y = tajima_D, fill = population)) +
    geom_col(alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    scale_x_discrete(labels = site_labels_stacked) +
    labs(x = NULL, y = expression(paste("Tajima's ", italic(D)))) +
    theme_report(11) + theme(legend.position = "none")
ggsave(file.path(plot_dir, "05_tajimasD_by_population.pdf"), p, width = 7, height = 6)
cat("Tajima's D plot written\n")
REOF

# -----------------------------------------------------------------------
# Export to PNG (300dpi) and SVG. pdftoppm/pdf2svg live in owenspopgen.
# -----------------------------------------------------------------------
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
command -v pdftoppm >/dev/null 2>&1 || { echo "ERROR: pdftoppm not found in owenspopgen env"; exit 1; }
command -v pdf2svg  >/dev/null 2>&1 || { echo "ERROR: pdf2svg not found in owenspopgen env"; exit 1; }

pdf="$PLOT_DIR/05_tajimasD_by_population.pdf"
if [[ -f "$pdf" ]]; then
    base="${pdf%.pdf}"
    pdftoppm -png -r 300 -singlefile "$pdf" "$base"
    pdf2svg "$pdf" "${base}.svg"
    echo "  exported: ${base}.png ${base}.svg"
fi
