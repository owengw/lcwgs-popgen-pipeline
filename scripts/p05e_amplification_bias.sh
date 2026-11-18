#!/bin/bash
#SBATCH --job-name=amp_bias_glm
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/amp_bias_glm_%j.log

# =============================================================================
# p05e_amplification_bias.sh — paired binomial GLM testing for WGA
# amplification bias in genome coverage breadth, using the 10 confirmed
# amplified-replicate / original-sample pairs (matched programmatically
# and verified against physalia/depths/all_depth_stats.tsv).
#
# Response: bases covered vs. genome size (binomial), i.e. breadth of
# coverage. Individual included as a blocking factor to respect the
# paired design (each pair shares one true individual).
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen   # R lives in this env
module load R/4.4.1-foss-2022b 2>/dev/null || true
set -eu

DEPTH_FILE="physalia/depths/all_depth_stats.tsv"
OUT_DIR="physalia/angsd_results/amp_bias"
mkdir -p "$OUT_DIR"

Rscript --vanilla - "$DEPTH_FILE" "$OUT_DIR" << 'REOF'
args       <- commandArgs(trailingOnly = TRUE)
depth_file <- args[1]
out_dir    <- args[2]

dat <- read.table(depth_file, header = FALSE, stringsAsFactors = FALSE)
colnames(dat)[1:8] <- c("sample_full","genome_size","mean_depth","col4",
                        "col5","bases_covered","breadth","col8")
dat$genome_size   <- as.numeric(dat$genome_size)
dat$bases_covered <- as.numeric(dat$bases_covered)
dat$breadth       <- as.numeric(dat$breadth)
dat$mean_depth    <- as.numeric(dat$mean_depth)

get_row <- function(pattern) {
    hit <- dat[grepl(pattern, dat$sample_full), ]
    if (nrow(hit) != 1) stop("Expected 1 match for ", pattern, ", got ", nrow(hit))
    hit
}

# confirmed pairs: (original unamplified DNA "o", amplified library of the same extract, individual_id)
# The T / aT re-extractions are not used: they are amplified too, so they cannot test amplification.
pairs <- list(
    c("289o-NS40o_merged", "283-NS40_merged",  "NS40"),
    c("294o-SI7o_merged",  "221-SI7_merged",   "SI7"),
    c("296o-NS14o_merged", "227-NS14_merged",  "NS14"),
    c("297o-NS41o_merged", "291-NS41_merged",  "NS41"),
    c("315o-SI18o_merged", "198-SI18_merged",  "SI18"),
    c("332o-CAI6o_merged", "327-CAI6_merged",  "CAI6"),
    c("337o-SI59o_merged", "233-SI59_merged",  "SI59"),
    c("351o-NS1o_merged",  "258-NS1_merged",   "NS1"),
    c("352o-NS30o_merged", "236-NS30_merged",  "NS30"),
    c("359o-NS32o_merged", "244-NS32_merged",  "NS32")
)

rows <- do.call(rbind, lapply(pairs, function(p) {
    raw_row <- get_row(p[1]); amp_row <- get_row(p[2])
    rbind(
        data.frame(individual = p[3], amplified = "amplified",
                   bases_covered = amp_row$bases_covered,
                   genome_size = amp_row$genome_size,
                   breadth = amp_row$breadth, mean_depth = amp_row$mean_depth),
        data.frame(individual = p[3], amplified = "non_amplified",
                   bases_covered = raw_row$bases_covered,
                   genome_size = raw_row$genome_size,
                   breadth = raw_row$breadth, mean_depth = raw_row$mean_depth)
    )
}))

cat("=== Paired coverage data (n=10 individuals) ===\n")
print(rows[, c("individual","amplified","breadth","mean_depth")])

rows$amplified <- factor(rows$amplified, levels = c("non_amplified", "amplified"))
rows$individual <- factor(rows$individual)
rows$not_covered <- rows$genome_size - rows$bases_covered

cat("\n=== Paired binomial GLM: coverage breadth ~ amplification status + individual ===\n")
mod <- glm(cbind(bases_covered, not_covered) ~ amplified + individual,
           data = rows, family = binomial)
print(summary(mod)$coefficients["amplifiedamplified", ])

cat("\n=== Mean breadth by group ===\n")
print(tapply(rows$breadth, rows$amplified, mean))

cat("\n=== Paired difference (amplified - non_amplified), per individual ===\n")
wide <- reshape(rows[, c("individual","amplified","breadth")],
                idvar = "individual", timevar = "amplified", direction = "wide")
wide$diff <- wide$breadth.amplified - wide$breadth.non_amplified
print(wide[, c("individual", "breadth.non_amplified", "breadth.amplified", "diff")])
cat("\nMean paired difference:", round(mean(wide$diff), 4), "\n")
cat("Wilcoxon signed-rank test (non-parametric check):\n")
print(wilcox.test(wide$breadth.amplified, wide$breadth.non_amplified, paired = TRUE))

write.table(rows, file.path(out_dir, "amplification_bias_data.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(mod, file.path(out_dir, "amplification_bias_glm.rds"))
cat("\nData and model written to:", out_dir, "\n")
REOF
