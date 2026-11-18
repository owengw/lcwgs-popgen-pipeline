#!/bin/bash
#SBATCH --job-name=gone_summary
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/gone_summary_%j.log

# =============================================================================
# p_gone_summary.sh — reads GONE's Output_Ne_<POP> files (CAI/NS/SB/SI) and
# produces (1) a summary table of Ne per population and (2) a line plot of
# Ne vs generation, one line per population.
#
# GONE FORMAT NOTE: Output_Ne_<POP> has a title line ("Ne averages over N
# independent estimates."), then a header ("Generation\tGeometric_mean"),
# then the table. The first few generations (here: 1-4) are commonly flat
# at an identical value - this is a known GONE artifact: LD signal isn't
# resolved into separate bins at the very shortest recombination distances,
# so those generations collapse to one estimate rather than reflecting
# genuinely constant Ne. The table below flags this explicitly and reports
# a "recent Ne" window that starts AFTER the flat lead-in, rather than
# quoting generation 1 as if it were a clean contemporary estimate.
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -eu

NE_DIR=${NE_DIR:-/mnt/parscratch/users/bi4og/physalia/ne_estimate}
PLOT_DIR="$NE_DIR/gone_summary"
mkdir -p "$PLOT_DIR"

Rscript --vanilla - "$NE_DIR" "$PLOT_DIR" << 'REOF'
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

theme_report <- function(base_size = 13) {
  theme_classic(base_size = base_size) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.key = element_rect(fill = "white", colour = NA))
}

args <- commandArgs(trailingOnly = TRUE)
ne_dir <- args[1]; plot_dir <- args[2]

site_labels_full <- c(CAI = "Castle Island", NS = "Nonsuch Island",
                      SB = "Sinky Bay",       SI = "Southampton Island")
pops <- names(site_labels_full)

# ---- read each Output_Ne_<POP> file (skip the title line, keep the header) ----
read_gone <- function(pop) {
  f <- file.path(ne_dir, pop, paste0("Output_Ne_", pop))
  if (!file.exists(f)) { cat("WARNING:", f, "not found - skipping", pop, "\n"); return(NULL) }
  dt <- fread(f, skip = 1, header = TRUE)  # skip "Ne averages over N..." title line
  setnames(dt, c("Generation", "Ne"))
  dt[, population := pop]
  dt
}
ne_all <- rbindlist(lapply(pops, read_gone))
if (nrow(ne_all) == 0) stop("No GONE output files found - check NE_DIR")

fwrite(ne_all, file.path(plot_dir, "gone_Ne_all_generations.tsv"), sep = "\t")

# ---- detect the flat lead-in per population (consecutive identical Ne from gen 1) ----
detect_flat_leadin <- function(dt) {
  dt <- dt[order(Generation)]
  first_val <- dt$Ne[1]
  run <- which(dt$Ne != first_val)[1]
  if (is.na(run)) length(dt$Ne) else run - 1L   # number of flat generations at the start
}
flat_len <- ne_all[, .(flat_leadin_generations = detect_flat_leadin(.SD)), by = population]

# ---- summary table: Ne at gen 1 (flagged), and "recent Ne" = mean over the
# first RECENT_WINDOW generations AFTER the flat lead-in ends ----
RECENT_WINDOW <- 5
summary_tab <- merge(
  ne_all[Generation == 1, .(population, Ne_generation1 = Ne)],
  flat_len, by = "population")

summary_tab[, recent_window_start := flat_leadin_generations + 1L]
summary_tab[, recent_window_end   := flat_leadin_generations + RECENT_WINDOW]

recent_ne <- ne_all[summary_tab, on = "population",
                     .(population, x.Generation, x.Ne,
                       recent_window_start = i.recent_window_start,
                       recent_window_end   = i.recent_window_end)][
    x.Generation >= recent_window_start & x.Generation <= recent_window_end,
    .(Ne_recent_mean = mean(x.Ne)), by = population]

summary_tab <- merge(summary_tab, recent_ne, by = "population")
summary_tab[, site := site_labels_full[population]]
summary_tab <- summary_tab[, .(site, population, flat_leadin_generations,
                                Ne_generation1, recent_window_start,
                                recent_window_end, Ne_recent_mean)]

cat("=== GONE Ne summary ===\n")
cat("NOTE: 'flat_leadin_generations' = how many generations from gen 1 share an\n")
cat("identical Ne estimate (a known GONE artifact, not genuine constancy).\n")
cat("Ne_recent_mean is the mean Ne over the", RECENT_WINDOW, "generations\n")
cat("immediately AFTER that flat lead-in - a more defensible 'contemporary Ne'\n")
cat("than quoting generation 1 directly.\n\n")
print(summary_tab)
fwrite(summary_tab, file.path(plot_dir, "gone_Ne_summary_table.tsv"), sep = "\t")

# ---- line plot: Ne vs generation, one line per population ----
ne_all[, site := site_labels_full[population]]
p <- ggplot(ne_all, aes(x = Generation, y = Ne, colour = site)) +
    geom_line(linewidth = 0.9) +
    labs(x = "Generations before present", y = expression(italic(N)[e]), colour = "Site") +
    theme_report(13)
ggsave(file.path(plot_dir, "gone_Ne_by_generation.pdf"), p, width = 8, height = 5.5)

p_log <- p + scale_y_log10()
ggsave(file.path(plot_dir, "gone_Ne_by_generation_logscale.pdf"), p_log, width = 8, height = 5.5)

cat("\nPlots and table written to:", plot_dir, "\n")
REOF

# -----------------------------------------------------------------------
# Export PDFs to PNG (300dpi) and SVG.
# -----------------------------------------------------------------------
command -v pdftoppm >/dev/null 2>&1 || { echo "ERROR: pdftoppm not found in owenspopgen env"; exit 1; }
command -v pdf2svg  >/dev/null 2>&1 || { echo "ERROR: pdf2svg not found in owenspopgen env"; exit 1; }

for pdf in "$PLOT_DIR"/*.pdf; do
    [[ -f "$pdf" ]] || continue
    base="${pdf%.pdf}"
    pdftoppm -png -r 300 -singlefile "$pdf" "$base"
    pdf2svg "$pdf" "${base}.svg"
    echo "  exported: ${base}.png ${base}.svg"
done

echo ""
echo "Done. Table + plots in: $PLOT_DIR"
