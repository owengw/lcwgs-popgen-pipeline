#!/bin/bash
#SBATCH --job-name=ngsrelate
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p07c_ngsrelate_%j.log

# =============================================================================
# ngsRelate pairwise relatedness estimation
# Fixed version:
#   - Skips beagle generation if file already exists
#   - Skips ngsRelate if .res file already exists
#   - Correct column names (theta not Kin, a/b not ida/idb)
#   - Appends to summary file rather than overwriting
#   - Clean sample ID extraction from BAM paths
#   - Base R only — no tidyverse dependency
#   - module load R for cluster compatibility
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
module load R/4.4.1-foss-2022b 2>/dev/null || true

set -eu
set -x

# =============================================================================
# Variables
# =============================================================================
BAM_LIST=${BAM_LIST:-}
REFERENCE=${REFERENCE:-}
OUT_DIR=${OUT_DIR:-}
THREADS=${THREADS:-8}

echo "=================================================================="
echo "ngsRelate: Pairwise Relatedness Estimation"
echo "=================================================================="
echo "BAM list:         $BAM_LIST"
echo "Output directory: $OUT_DIR"
echo "=================================================================="

mkdir -p "${OUT_DIR}/relatedness"

POP_NAME=$(basename "$BAM_LIST" .list)
OUT_PREFIX="${OUT_DIR}/relatedness/relate_${POP_NAME}"
BEAGLE_FILE="${OUT_PREFIX}_relate.beagle.gz"
RES_FILE="${OUT_PREFIX}.res"

N_SAMPLES=$(wc -l < "$BAM_LIST")
echo "Population: $POP_NAME"
echo "Sample size: $N_SAMPLES"

# Skip small populations — ngsRelate is not meaningful below n=5
if [[ $N_SAMPLES -lt 5 ]]; then
    echo "WARNING: $POP_NAME has only $N_SAMPLES samples — skipping ngsRelate"
    exit 0
fi

# =============================================================================
# Step 1 — Generate beagle file (skip if already exists)
# =============================================================================

if [[ -f "$BEAGLE_FILE" ]]; then
    echo "Beagle file already exists, skipping ANGSD: $BEAGLE_FILE"
else
    echo "Generating beagle file from random 10M sites..."

    SITES_FILE="${OUT_DIR}/random_10M_sites.list"
    RF_FILE="${OUT_DIR}/random_10M_contigs.txt"

    if [[ ! -f "$SITES_FILE" ]]; then
        echo "ERROR: Random sites file not found: $SITES_FILE"
        exit 1
    fi

    MIN_IND=$(echo "$N_SAMPLES * 0.2" | bc | awk '{print int($1+0.5)}')
    [[ $MIN_IND -lt 1 ]] && MIN_IND=1
    echo "minInd: $MIN_IND"

    angsd \
        -bam "$BAM_LIST" \
        -ref "$REFERENCE" \
        -GL 2 \
        -doGlf 2 \
        -doMajorMinor 1 \
        -doMaf 1 \
        -minMapQ 10 \
        -minQ 20 \
        -minInd $MIN_IND \
        -sites "$SITES_FILE" \
        -rf "$RF_FILE" \
        -P $THREADS \
        -out "${OUT_PREFIX}_relate" \
        2> "${OUT_PREFIX}_beagle.log"

    if [[ ! -f "$BEAGLE_FILE" ]]; then
        echo "ERROR: Failed to generate beagle file"
        exit 1
    fi

    echo "Beagle file generated: $BEAGLE_FILE"
fi

# Get number of sites from beagle file
N_SITES=$(zcat "$BEAGLE_FILE" | tail -n +2 | wc -l)
echo "Sites in beagle: $N_SITES"

# =============================================================================
# Step 2 — Run ngsRelate (skip if .res already exists)
# =============================================================================

if [[ -f "$RES_FILE" ]]; then
    echo "ngsRelate output already exists, skipping: $RES_FILE"
else
    echo "Running ngsRelate..."

    ngsRelate \
        -G "$BEAGLE_FILE" \
        -n $N_SAMPLES \
        -L $N_SITES \
        -O "$RES_FILE" \
        -p $THREADS \
        2> "${OUT_PREFIX}.log"

    if [[ ! -f "$RES_FILE" ]]; then
        echo "ERROR: ngsRelate failed — check ${OUT_PREFIX}.log"
        tail -50 "${OUT_PREFIX}.log"
        exit 1
    fi

    echo "ngsRelate complete: $RES_FILE"
fi

# =============================================================================
# Step 3 — Process results and generate summary
# Base R only — no tidyverse dependency
# =============================================================================

echo "Processing relatedness results..."

Rscript --vanilla - "$RES_FILE" "$BAM_LIST" "$OUT_PREFIX" << 'REOF'

args       <- commandArgs(trailingOnly = TRUE)
res_file   <- args[1]
bam_list   <- args[2]
out_prefix <- args[3]

# -------------------------------------------------------
# Load results
# -------------------------------------------------------
relate     <- read.table(res_file, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
bam_paths  <- readLines(bam_list)
sample_ids <- sub("_merged.*", "", basename(bam_paths))

cat("Samples loaded:", length(sample_ids), "\n")
cat("Pairs in .res file:", nrow(relate), "\n")

# ngsRelate uses 0-based indexing for a and b columns
relate$sample_a <- sample_ids[relate$a + 1]
relate$sample_b <- sample_ids[relate$b + 1]

# Reorder so sample_a and sample_b are first columns
other_cols   <- setdiff(names(relate), c("sample_a", "sample_b"))
relate       <- relate[, c("sample_a", "sample_b", other_cols)]

# Write full results with sample names
full_file <- paste0(out_prefix, "_full.txt")
write.table(relate, full_file,
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Full results written to:", full_file, "\n")

# -------------------------------------------------------
# Classify relationships using KING coefficient
# KING thresholds (Manichaikul et al. 2010):
#   > 0.354  : duplicate / monozygotic twin
#   > 0.177  : 1st degree (parent-offspring or full siblings)
#   > 0.0884 : 2nd degree (half-siblings, grandparent-grandchild)
#   > 0.0442 : 3rd degree (first cousins)
#   <= 0.0442: unrelated
# -------------------------------------------------------
close_rels <- relate[relate$KING > 0.0442, ]
close_rels <- close_rels[order(-close_rels$KING), ]

close_rels$relationship <- ifelse(
    close_rels$KING > 0.354,  "Duplicate/Identical twin",
    ifelse(close_rels$KING > 0.177,  "1st degree (parent-offspring/full sib)",
    ifelse(close_rels$KING > 0.0884, "2nd degree (half-sib/grandparent)",
                                     "3rd degree (first cousin)")))

close_file <- paste0(out_prefix, "_close_relationships.txt")
write.table(close_rels, close_file,
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Close relationships written to:", close_file, "\n")

# -------------------------------------------------------
# Print close relationships to console
# -------------------------------------------------------
pop_name  <- basename(out_prefix)
n_samples <- length(sample_ids)
n_pairs   <- nrow(relate)

cat("\n=== Close relationships for", pop_name, "===\n")
if (nrow(close_rels) > 0) {
    display <- close_rels[, c("sample_a", "sample_b", "nSites",
                               "KING", "relationship")]
    display$sample_a <- sub("^[0-9]+-", "", display$sample_a)
    display$sample_b <- sub("^[0-9]+-", "", display$sample_b)
    display$KING     <- round(display$KING, 4)
    print(display, row.names = FALSE)
} else {
    cat("  None detected above KING threshold of 0.0442\n")
}

# -------------------------------------------------------
# Append to shared summary file — one block per population
# -------------------------------------------------------
summary_file <- file.path(dirname(out_prefix), "relatedness_summary.txt")

if (!file.exists(summary_file)) {
    writeLines(c("Relatedness Analysis Summary",
                 "============================="),
               summary_file)
}

lines <- c(
    "",
    paste0("Population: ", pop_name),
    paste0("  Number of samples:               ", n_samples),
    paste0("  Number of pairs:                 ", n_pairs),
    paste0("  Mean KING coefficient:           ",
           round(mean(relate$KING,  na.rm = TRUE), 6)),
    paste0("  Mean theta (kinship):            ",
           round(mean(relate$theta, na.rm = TRUE), 6)),
    "  Close relationships:",
    paste0("    Duplicates/identical twins (KING > 0.354):    ",
           sum(close_rels$KING > 0.354)),
    paste0("    1st degree (KING 0.177-0.354):                ",
           sum(close_rels$KING > 0.177  & close_rels$KING <= 0.354)),
    paste0("    2nd degree (KING 0.0884-0.177):               ",
           sum(close_rels$KING > 0.0884 & close_rels$KING <= 0.177)),
    paste0("    3rd degree (KING 0.0442-0.0884):              ",
           sum(close_rels$KING > 0.0442 & close_rels$KING <= 0.0884))
)

cat(paste(lines, collapse = "\n"), "\n",
    file = summary_file, append = TRUE)
cat("\nSummary appended to:", summary_file, "\n")
cat("Complete!\n")
REOF

echo ""
echo "=================================================================="
echo "ngsRelate complete: $POP_NAME"
echo "  Full results:        ${OUT_PREFIX}_full.txt"
echo "  Close relationships: ${OUT_PREFIX}_close_relationships.txt"
echo "  Summary:             ${OUT_DIR}/relatedness/relatedness_summary.txt"
echo "=================================================================="