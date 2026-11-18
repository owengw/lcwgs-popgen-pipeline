#!/bin/bash
#SBATCH --job-name=collate_downsampled
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/p07d_collate_downsampled_%j.log

# =============================================================================
# p07d_collate_downsampled.sh — run after all p07d_theta_downsampled.sh chains
# complete. Compares downsampled (n=16) vs full-n Tajima's D per population.
# Base R only.
# =============================================================================

source ~/.bash_profile
module load R/4.4.1-foss-2022b 2>/dev/null || true
set -eu

OUT_DIR=${OUT_DIR:-physalia/angsd_results}
DS_DIR="$OUT_DIR/theta_downsampled"
THETA_DIR="$OUT_DIR/theta_corrected"
DOWNSAMPLE_N=16
OUT_FILE="$DS_DIR/theta_downsampled_summary.tsv"

Rscript --vanilla - "$OUT_DIR" "$THETA_DIR" "$DOWNSAMPLE_N" "$OUT_FILE" << 'REOF'
args     <- commandArgs(trailingOnly = TRUE)
out_dir  <- args[1]
full_dir <- args[2]
ds_n     <- as.integer(args[3])
out_file <- args[4]

populations <- c("CAI", "NS", "SB", "SI")

read_pestpg <- function(path, pop, label) {
    if (!file.exists(path)) { cat("WARNING: missing", path, "\n"); return(NULL) }
    dat <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
    if (nrow(dat) == 0) return(NULL)
    n_sites  <- sum(dat$nSites, na.rm = TRUE)
    theta_pi <- sum(dat$tP,     na.rm = TRUE)
    theta_w  <- sum(dat$tW,     na.rm = TRUE)
    tajima_D <- weighted.mean(dat$Tajima, dat$nSites, na.rm = TRUE)
    data.frame(population = pop, analysis = label, n_sites = n_sites,
               theta_pi = theta_pi, theta_w = theta_w, tajima_D = tajima_D,
               theta_pi_per_site = theta_pi / n_sites,
               theta_w_per_site  = theta_w  / n_sites,
               stringsAsFactors = FALSE)
}

ds_summary <- do.call(rbind, lapply(populations, function(p)
    read_pestpg(file.path(out_dir, sprintf("pop_map_%s_ds%d.pestPG", p, ds_n)),
                p, sprintf("downsampled_n%d", ds_n))))

full_summary <- do.call(rbind, lapply(populations, function(p)
    read_pestpg(file.path(out_dir, sprintf("pop_map_%s.pestPG", p)),
                p, "full_n")))

combined <- rbind(ds_summary, full_summary)
combined <- combined[order(combined$population, combined$analysis), ]
write.table(combined, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Downsampled vs full-n theta comparison ===\n")
print(combined[, c("population", "analysis", "n_sites",
                   "theta_pi_per_site", "theta_w_per_site", "tajima_D")])

cat("\n=== Tajima's D comparison (key output) ===\n")
tajima_wide <- reshape(combined[, c("population", "analysis", "tajima_D")],
                       idvar = "population", timevar = "analysis", direction = "wide")
names(tajima_wide) <- sub("^tajima_D\\.", "", names(tajima_wide))
print(tajima_wide)

cat("\nResults written to:", out_file, "\n")
REOF
