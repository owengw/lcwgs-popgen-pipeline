#!/bin/bash
#SBATCH --job-name=theta_downsampled
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p07d_theta_downsampled_%j.log

# =============================================================================
# Downsampled theta analysis — equal n=17 per population
# Populations: CAI, NS, SB, SI
# Purpose: remove sample size bias from Tajima's D comparison
#
# Fixed version:
#   - stdout and stderr separated for realSFS (tee was corrupting .sfs files)
#   - SAF skip logic to avoid rerunning slow ANGSD if SAF already exists
#   - SFS entry count validation before proceeding to saf2theta
#   - Subsampling skip logic if lists already exist
#   - R collation uses base R only (no tidyverse dependency)
#   - module load R for cluster compatibility
#
# SB already has exactly 17 clean individuals so is used as-is
# SI83 excluded as duplicate of SI77
# Random seed = 42 — report in methods for reproducibility
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
module load R/4.4.1-foss-2022b 2>/dev/null || true

set -euo pipefail
set -x

# =============================================================================
# Variables
# =============================================================================
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
REFERENCE="genome/Plong_genome_flye.fasta"
SITES_FILE="$OUT_DIR/random_10M_sites.list"
RF_FILE="$OUT_DIR/random_10M_contigs.txt"
POP_LIST_DIR="$OUT_DIR/pop_map_lists"
THETA_DIR="$OUT_DIR/theta_downsampled"
THREADS=6
DOWNSAMPLE_N=17
SEED=42
WINDOW=50000
STEP=25000
MIN_IND=5   # 30% of n=17

mkdir -p "$THETA_DIR"
mkdir -p "$THETA_DIR/subsampled_lists"

POPULATIONS=(CAI NS SB SI)

# Expected SFS entries for n=17 folded: 2*17+1 = 35
EXPECTED_SFS_ENTRIES=35

# =============================================================================
# Step 1 — Create clean subsampled BAM lists
# Skip if subsampled list already exists
# =============================================================================

echo "Creating subsampled BAM lists (n=$DOWNSAMPLE_N, seed=$SEED)..."

for POP in "${POPULATIONS[@]}"; do
    FULL_LIST="$POP_LIST_DIR/pop_map_${POP}.list"
    CLEAN_LIST="$THETA_DIR/subsampled_lists/${POP}_clean.list"
    SUB_LIST="$THETA_DIR/subsampled_lists/${POP}_n${DOWNSAMPLE_N}.list"

    if [[ ! -f "$FULL_LIST" ]]; then
        echo "ERROR: BAM list not found: $FULL_LIST"
        exit 1
    fi

    if [[ -f "$SUB_LIST" ]]; then
        echo "  $POP: subsampled list already exists, skipping: $SUB_LIST"
        continue
    fi

    # Remove duplicates, originals, and the SI83 known duplicate
    grep -v "o_merged\|T_merged\|SI83_merged" "$FULL_LIST" > "$CLEAN_LIST"

    N_CLEAN=$(wc -l < "$CLEAN_LIST")
    echo "  $POP: $N_CLEAN clean individuals"

    if [[ $N_CLEAN -lt $DOWNSAMPLE_N ]]; then
        echo "ERROR: $POP has only $N_CLEAN clean individuals, need $DOWNSAMPLE_N"
        exit 1
    fi

    if [[ $N_CLEAN -eq $DOWNSAMPLE_N ]]; then
        cp "$CLEAN_LIST" "$SUB_LIST"
        echo "  $POP: already at n=$DOWNSAMPLE_N, using full clean list"
    else
        shuf --random-source=<(openssl enc -aes-256-ctr \
            -pass pass:"seed${SEED}" -nosalt \
            </dev/zero 2>/dev/null) \
            -n $DOWNSAMPLE_N "$CLEAN_LIST" > "$SUB_LIST"
        echo "  $POP: subsampled $N_CLEAN -> $DOWNSAMPLE_N individuals"
    fi

    echo "  Subsampled list written: $SUB_LIST"
    echo ""
done

# =============================================================================
# Step 2 — Run theta pipeline for each subsampled population
# =============================================================================

for POP in "${POPULATIONS[@]}"; do
    echo "=========================================="
    echo "Processing population: $POP (n=$DOWNSAMPLE_N)"
    echo "=========================================="

    BAM_LIST="$THETA_DIR/subsampled_lists/${POP}_n${DOWNSAMPLE_N}.list"
    OUT_PREFIX="$THETA_DIR/pop_map_${POP}_n${DOWNSAMPLE_N}"

    # ----------------------------------------------------------
    # SAF generation — skip if already exists
    # ----------------------------------------------------------
    if [[ -f "${OUT_PREFIX}.saf.idx" ]]; then
        echo "  SAF already exists, skipping ANGSD: ${OUT_PREFIX}.saf.idx"
    else
        angsd \
            -bam "$BAM_LIST" \
            -ref "$REFERENCE" \
            -anc "$REFERENCE" \
            -GL 2 \
            -doSaf 1 \
            -doCounts 1 \
            -minMapQ 10 \
            -minQ 20 \
            -minInd $MIN_IND \
            -setMinDepthInd 1 \
            -P $THREADS \
            -nQueueSize 50 \
            -remove_bads 1 \
            -trim 0 \
            -sites "$SITES_FILE" \
            -rf "$RF_FILE" \
            -out "$OUT_PREFIX" \
            2> "${OUT_PREFIX}.angsd_saf.log"
    fi

    # ----------------------------------------------------------
    # Folded SFS
    # CRITICAL: stdout to .sfs ONLY, stderr to .sfs.log separately
    # Do NOT use 2>&1 | tee — this pipes log messages into the SFS
    # file and corrupts it, causing saf2theta to fail with dimension errors
    # Validate entry count before proceeding
    # ----------------------------------------------------------
    ACTUAL_SFS_ENTRIES=0
    [[ -f "${OUT_PREFIX}.sfs" ]] && ACTUAL_SFS_ENTRIES=$(wc -w < "${OUT_PREFIX}.sfs")

    if [[ $ACTUAL_SFS_ENTRIES -eq $EXPECTED_SFS_ENTRIES ]]; then
        echo "  Valid SFS already exists ($ACTUAL_SFS_ENTRIES entries), skipping realSFS"
    else
        echo "  Running realSFS (expected $EXPECTED_SFS_ENTRIES entries)..."
        realSFS \
            "${OUT_PREFIX}.saf.idx" \
            -fold 1 \
            -P $THREADS \
            > "${OUT_PREFIX}.sfs" \
            2> "${OUT_PREFIX}.sfs.log"

        ACTUAL_SFS_ENTRIES=$(wc -w < "${OUT_PREFIX}.sfs")
        echo "  SFS entries produced: $ACTUAL_SFS_ENTRIES (expected $EXPECTED_SFS_ENTRIES)"

        if [[ $ACTUAL_SFS_ENTRIES -ne $EXPECTED_SFS_ENTRIES ]]; then
            echo "ERROR: SFS has wrong number of entries for $POP"
            echo "  Expected: $EXPECTED_SFS_ENTRIES  Got: $ACTUAL_SFS_ENTRIES"
            echo "  Check ${OUT_PREFIX}.sfs.log for errors"
            exit 1
        fi
    fi

    # ----------------------------------------------------------
    # Theta estimation from SFS
    # stderr to log, stdout goes to terminal via set -x
    # ----------------------------------------------------------
    realSFS saf2theta \
        "${OUT_PREFIX}.saf.idx" \
        -sfs "${OUT_PREFIX}.sfs" \
        -fold 1 \
        -outname "$OUT_PREFIX" \
        2> "${OUT_PREFIX}.saf2theta.log"

    # ----------------------------------------------------------
    # Windowed theta statistics
    # ----------------------------------------------------------
    thetaStat do_stat \
        "${OUT_PREFIX}.thetas.idx" \
        -win $WINDOW \
        -step $STEP \
        -outnames "${OUT_PREFIX}.thetas.windowed" \
        2> "${OUT_PREFIX}.thetas.windowed.log"

    # ----------------------------------------------------------
    # Genome-wide summary
    # stdout to .thetas.pestPG ONLY, stderr to log
    # ----------------------------------------------------------
    thetaStat print \
        "${OUT_PREFIX}.thetas.idx" \
        > "${OUT_PREFIX}.thetas.pestPG" \
        2> "${OUT_PREFIX}.thetas.stats.log"

    PESTPG_LINES=$(wc -l < "${OUT_PREFIX}.thetas.pestPG")
    echo "  pestPG lines: $PESTPG_LINES"
    echo "  Done: $POP"
done

# =============================================================================
# Step 3 — Collate and compare with full-n results
# Base R only — no tidyverse dependency
# =============================================================================

echo "Collating downsampled theta summary..."

Rscript --vanilla - \
    "$THETA_DIR" \
    "$OUT_DIR/theta_corrected" \
    "$DOWNSAMPLE_N" \
    "$THETA_DIR/theta_downsampled_summary.tsv" << 'REOF'

args     <- commandArgs(trailingOnly = TRUE)
ds_dir   <- args[1]
full_dir <- args[2]
ds_n     <- as.integer(args[3])
out_file <- args[4]

# Use windowed pestPG — has tW, tP, Tajima, nSites columns
# comment.char="" required — header line starts with #
# Per-site pestPG (thetaStat print) has different columns and is not used
read_pestpg <- function(dir, pattern, label) {
    files <- list.files(dir, pattern = pattern, full.names = TRUE)
    files <- files[!grepl("pop_map_DH", files)]
    if (length(files) == 0) {
        cat("WARNING: No files matching", pattern, "in", dir, "\n")
        return(NULL)
    }
    results <- lapply(files, function(f) {
        pop <- sub("\\.thetas\\.windowed\\.pestPG$", "", basename(f))
        pop <- sub("_n[0-9]+$", "", pop)
        pop <- sub("^pop_map_", "", pop)

        dat <- tryCatch(
            read.table(f, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE,
                       comment.char = ""),
            error = function(e) {
                cat("ERROR reading", f, ":", conditionMessage(e), "\n")
                return(NULL)
            }
        )
        if (is.null(dat) || nrow(dat) == 0) {
            cat("WARNING: Empty file:", f, "\n")
            return(NULL)
        }

        n_sites  <- sum(dat$nSites, na.rm = TRUE)
        theta_pi <- sum(dat$tP,     na.rm = TRUE)
        theta_w  <- sum(dat$tW,     na.rm = TRUE)
        tajima_D <- weighted.mean(dat$Tajima, dat$nSites, na.rm = TRUE)

        data.frame(
            population        = pop,
            analysis          = label,
            n_sites           = n_sites,
            theta_pi          = theta_pi,
            theta_w           = theta_w,
            tajima_D          = tajima_D,
            theta_pi_per_site = theta_pi / n_sites,
            theta_w_per_site  = theta_w  / n_sites,
            stringsAsFactors  = FALSE
        )
    })
    results <- results[!sapply(results, is.null)]
    if (length(results) == 0) return(NULL)
    do.call(rbind, results)
}

ds_summary <- read_pestpg(
    ds_dir,
    pattern = "\\.thetas\\.windowed\\.pestPG$",
    label   = paste0("downsampled_n", ds_n)
)

full_summary <- read_pestpg(
    full_dir,
    pattern = "\\.thetas\\.windowed\\.pestPG$",
    label   = "full_n"
)
if (!is.null(full_summary)) {
    full_summary <- full_summary[
        full_summary$population %in% c("CAI", "NS", "SB", "SI"), ]
}

combined <- rbind(ds_summary, full_summary)
combined  <- combined[order(combined$population, combined$analysis), ]

write.table(combined, out_file,
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Downsampled vs full-n theta comparison ===\n")
print(combined[, c("population", "analysis", "n_sites",
                   "theta_pi_per_site", "theta_w_per_site", "tajima_D")])

cat("\n=== Tajima's D comparison (key output) ===\n")
tajima_wide <- reshape(
    combined[, c("population", "analysis", "tajima_D")],
    idvar     = "population",
    timevar   = "analysis",
    direction = "wide"
)
names(tajima_wide) <- sub("^tajima_D\\.", "", names(tajima_wide))
print(tajima_wide)

cat("\nResults written to:", out_file, "\n")
REOF

echo ""
echo "=========================================="
echo "Downsampled theta analysis complete"
echo "Results in: $THETA_DIR"
echo "Summary: $THETA_DIR/theta_downsampled_summary.tsv"
echo ""
echo "Subsampled BAM lists saved to:"
for POP in "${POPULATIONS[@]}"; do
    echo "  $POP: $THETA_DIR/subsampled_lists/${POP}_n${DOWNSAMPLE_N}.list"
done
echo "Random seed: $SEED"
echo "=========================================="