#!/bin/bash
#SBATCH --job-name=rohan_collate
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/p08b_rohan_collate_%j.log

# =============================================================================
# ROHan population-level collation
# Run AFTER p08a_rohan_individual.sh completes for all individuals
#
# Reads individual .summary files from roh/individual/
# Groups by population using metadata.tsv
# Produces population-level ROH summary statistics
#
# No ANGSD or ROHan required — pure R collation
# =============================================================================

source ~/.bash_profile
module load R/4.4.1-foss-2022b 2>/dev/null || true

set -euo pipefail

OUT_DIR=${OUT_DIR:-physalia/angsd_results}
ROH_IND_DIR="$OUT_DIR/roh/individual"
ROH_POP_DIR="$OUT_DIR/roh/population"
METADATA="$OUT_DIR/metadata.tsv"

mkdir -p "$ROH_POP_DIR"

echo "=================================================================="
echo "ROHan population collation"
echo "Individual results from: $ROH_IND_DIR"
echo "=================================================================="

# Check how many individual summaries exist
N_SUMMARIES=$(ls "$ROH_IND_DIR"/*.summary.txt 2>/dev/null | wc -l)
echo "Individual summary files found: $N_SUMMARIES"

if [[ $N_SUMMARIES -eq 0 ]]; then
    echo "ERROR: No individual summary files found in $ROH_IND_DIR"
    echo "Run p08a_rohan_individual.sh first and wait for all array tasks to complete"
    exit 1
fi

cat > /tmp/p08b_collate.R << 'REOF'

args        <- commandArgs(trailingOnly = TRUE)
ind_dir     <- args[1]
pop_dir     <- args[2]
meta_file   <- args[3]

# Load metadata for population assignment
metadata <- read.table(meta_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
# Clean sample prefix for matching
metadata$sample_prefix <- sub("_.*", "", metadata$sample_id)

# Find all individual summary files
summary_files <- list.files(ind_dir, pattern = "[.]summary[.]txt$", full.names = TRUE)

# Parse ROHan summary file
# ROHan .summary format varies slightly by version
# Key fields: heterozygosity (theta), number of ROH, total ROH length, F_ROH
parse_rohan_summary <- function(f) {
    sample <- sub("\\.summary$", "", basename(f))
    lines  <- tryCatch(readLines(f), error = function(e) NULL)
    if (is.null(lines) || length(lines) == 0) {
        cat("WARNING: Empty or unreadable summary:", f, "\n")
        return(NULL)
    }

    # Helper: extract first numeric value after label on matching line
    # ROHan format: "Label text :  VALUE (lower,upper)"
    get_val <- function(pattern) {
        m <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
        if (length(m) == 0) return(NA_real_)
        # Extract first number after the colon
        tok <- sub(".*:[[:space:]]*([0-9eE.+-]+).*", "\\1", m[1])
        suppressWarnings(as.numeric(tok))
    }

    data.frame(
        sample_id    = sample,
        theta        = get_val("Genome-wide theta outside ROH"),
        n_roh        = get_val("Segments in ROH[[:space:]]*:"),
        pct_roh      = get_val("Segments in ROH[(]"),
        avg_roh_bp   = get_val("Avg\\. length of ROH"),
        stringsAsFactors = FALSE
    )
}

ind_list <- lapply(summary_files, parse_rohan_summary)
ind_list <- ind_list[!sapply(ind_list, is.null)]

if (length(ind_list) == 0) {
    cat("ERROR: No summary files could be parsed\n")
    quit(status = 1)
}

ind_df <- do.call(rbind, ind_list)

# Show first parsed file raw to help diagnose column matching
cat("\nExample parsed values (first 5 individuals):\n")
print(head(ind_df, 5))

# Check for parsing failures (all NA)
n_failed <- sum(is.na(ind_df$theta) & is.na(ind_df$n_roh) &
                is.na(ind_df$pct_roh) & is.na(ind_df$avg_roh_bp))
if (n_failed > 0) {
    cat("\nWARNING:", n_failed, "individuals had no parseable values\n")
    cat("Check raw summary file format:\n")
    cat(readLines(summary_files[1]), sep = "\n")
}

# Match to population via metadata
ind_df$sample_prefix <- sub("_.*", "", ind_df$sample_id)
# Also handle numeric prefix format e.g. "14-SI88" -> prefix "14-SI88"
# Try direct match first, then strip numeric prefix
ind_df <- merge(ind_df,
                metadata[, c("sample_prefix", "site", "population", "year")],
                by = "sample_prefix", all.x = TRUE)

# Flag unmatched
n_unmatched <- sum(is.na(ind_df$site))
if (n_unmatched > 0) {
    cat("\nWARNING:", n_unmatched, "individuals not matched to metadata\n")
    print(ind_df[is.na(ind_df$site), c("sample_id", "sample_prefix")])
}

# Convert total ROH to Mb
ind_df$avg_roh_kb <- ind_df$avg_roh_bp / 1000

# Write individual results
ind_file <- file.path(pop_dir, "individual_ROH_all.tsv")
write.table(ind_df[order(ind_df$site, ind_df$sample_id), ],
            ind_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nIndividual ROH results written to:", ind_file, "\n")

# Population-level summary — main populations only
main_pops <- c("CAI", "NS", "SB", "SI")
ind_main  <- ind_df[!is.na(ind_df$site) & ind_df$site %in% main_pops, ]

pop_summary <- do.call(rbind, lapply(main_pops, function(pop) {
    d <- ind_main[ind_main$site == pop, ]
    if (nrow(d) == 0) return(NULL)
    data.frame(
        population        = pop,
        n_individuals     = nrow(d),
        mean_theta        = round(mean(d$theta,        na.rm = TRUE), 6),
        sd_theta          = round(sd(d$theta,          na.rm = TRUE), 6),
        mean_n_roh        = round(mean(d$n_roh,        na.rm = TRUE), 2),
        sd_n_roh          = round(sd(d$n_roh,          na.rm = TRUE), 2),
        mean_avg_roh_kb   = round(mean(d$avg_roh_kb, na.rm = TRUE), 3),
        sd_avg_roh_kb     = round(sd(d$avg_roh_kb,   na.rm = TRUE), 3),
        mean_pct_roh      = round(mean(d$pct_roh,    na.rm = TRUE), 4),
        sd_pct_roh        = round(sd(d$pct_roh,      na.rm = TRUE), 4),
        stringsAsFactors  = FALSE
    )
}))

pop_file <- file.path(pop_dir, "population_ROH_summary.tsv")
write.table(pop_summary, pop_file,
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Population ROH summary ===\n")
print(pop_summary)
cat("\nWritten to:", pop_file, "\n")

# Also include COI, SP, DH if present (with caveat about small n)
all_pops  <- unique(ind_df$site[!is.na(ind_df$site)])
small_pops <- setdiff(all_pops, main_pops)
if (length(small_pops) > 0) {
    cat("\nSmall populations (interpret with caution):\n")
    ind_small <- ind_df[!is.na(ind_df$site) & ind_df$site %in% small_pops, ]
    print(ind_small[, c("sample_id", "site", "theta", "n_roh",
    print(ind_small[, c("sample_id", "site", "theta", "n_roh", "avg_roh_kb", "pct_roh")])
}
REOF
Rscript --vanilla /tmp/p08b_collate.R "$ROH_IND_DIR" "$ROH_POP_DIR" "$METADATA"

echo ""
echo "=========================================="
echo "ROHan collation complete"
echo "Results in: $ROH_POP_DIR"
echo "=========================================="