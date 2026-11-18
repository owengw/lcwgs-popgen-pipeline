#!/bin/bash
#SBATCH --job-name=pop_theta
#SBATCH --cpus-per-task=6
#SBATCH --mem=250G
#SBATCH --time=72:00:00
#SBATCH --output=logs/p07b_pop_theta_%j.log

# =============================================================================
# Population-level theta/pi/Tajima's D estimation
# Fixed version:
#   - stdout and stderr separated for realSFS (tee was corrupting .sfs files)
#   - SAF skip logic to avoid rerunning slow ANGSD if SAF already exists
#   - DH skipped (n=1, population SFS not estimable)
#   - R collation uses base R only (no tidyverse dependency)
#   - module load R for cluster compatibility
#   - POP_SINGLE: optional variable to run a single population only
#     Usage: sbatch --export=OUT_DIR=physalia/angsd_results,POP_SINGLE=SI p07b_pop_theta.sh
#
# Supported populations:
#   CAI, COI, NS, SB, SI, SP  — per-population analyses
#   ALL                        — species-wide analysis (all_samples.list, minInd=59)
#
# ALL uses a separate BAM list: $OUT_DIR/all_samples.list
# Output goes to the same theta_corrected directory as population runs
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
THETA_DIR="$OUT_DIR/theta_corrected"
POP_LIST_DIR="$OUT_DIR/pop_map_lists"
THREADS=6
WINDOW=50000
STEP=25000

mkdir -p "$THETA_DIR"

# Populations to run — DH excluded (n=1, population SFS not estimable)
# ALL = species-wide run using all_samples.list
POPULATIONS=(CAI COI NS SB SI SP)

# If POP_SINGLE is set via --export, override to run only that population
# e.g. sbatch --export=OUT_DIR=physalia/angsd_results,POP_SINGLE=ALL p07b_pop_theta.sh
if [[ -n "${POP_SINGLE:-}" ]]; then
    echo "POP_SINGLE=$POP_SINGLE — running single population only"
    POPULATIONS=("$POP_SINGLE")
fi

# minInd thresholds per population (30% of n)
# ALL: 198 individuals total, 30% = 59
declare -A MIN_IND
MIN_IND[CAI]=7
MIN_IND[COI]=1
MIN_IND[NS]=16
MIN_IND[SB]=6
MIN_IND[SI]=28
MIN_IND[SP]=1
MIN_IND[ALL]=59

# BAM list paths — ALL uses a separate list outside pop_map_lists/
declare -A BAM_LISTS
BAM_LISTS[CAI]="$POP_LIST_DIR/pop_map_CAI.list"
BAM_LISTS[COI]="$POP_LIST_DIR/pop_map_COI.list"
BAM_LISTS[NS]="$POP_LIST_DIR/pop_map_NS.list"
BAM_LISTS[SB]="$POP_LIST_DIR/pop_map_SB.list"
BAM_LISTS[SI]="$POP_LIST_DIR/pop_map_SI.list"
BAM_LISTS[SP]="$POP_LIST_DIR/pop_map_SP.list"
BAM_LISTS[ALL]="$OUT_DIR/all_samples.list"

for POP in "${POPULATIONS[@]}"; do
    echo "=========================================="
    echo "Processing population: $POP"
    echo "=========================================="

    BAM_LIST="${BAM_LISTS[$POP]}"
    OUT_PREFIX="$THETA_DIR/pop_map_${POP}"

    if [[ ! -f "$BAM_LIST" ]]; then
        echo "WARNING: BAM list not found for $POP: $BAM_LIST — skipping"
        continue
    fi

    N_SAMPLES=$(wc -l < "$BAM_LIST")
    MIND=${MIN_IND[$POP]}
    echo "  Samples: $N_SAMPLES, minInd: $MIND"

    # ----------------------------------------------------------
    # Step 1 — Generate SAF
    # Skip if SAF index already exists to avoid rerunning slow ANGSD
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
            -minInd $MIND \
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
    # Step 2 — Estimate folded SFS
    # CRITICAL: stdout goes to .sfs ONLY, stderr goes to .sfs.log separately
    # Using 2>&1 | tee corrupts the .sfs file with log messages
    # Skip if valid SFS already exists (check entry count matches 2n+1)
    # ----------------------------------------------------------
    EXPECTED_SFS_ENTRIES=$(( N_SAMPLES * 2 + 1 ))
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
    # Step 3 — Estimate thetas from SFS
    # ----------------------------------------------------------
    realSFS saf2theta \
        "${OUT_PREFIX}.saf.idx" \
        -sfs "${OUT_PREFIX}.sfs" \
        -fold 1 \
        -outname "$OUT_PREFIX" \
        2> "${OUT_PREFIX}.saf2theta.log"

    # ----------------------------------------------------------
    # Step 4 — Windowed theta statistics
    # ----------------------------------------------------------
    thetaStat do_stat \
        "${OUT_PREFIX}.thetas.idx" \
        -win $WINDOW \
        -step $STEP \
        -outnames "${OUT_PREFIX}.thetas.windowed" \
        2> "${OUT_PREFIX}.thetas.windowed.log"

    # ----------------------------------------------------------
    # Step 5 — Genome-wide summary (per-site pestPG)
    # stdout goes to .thetas.pestPG ONLY, stderr to log
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
# Collate genome-wide theta summary — base R only, no tidyverse dependency
# Use windowed pestPG (has tW, tP, Tajima, nSites columns)
# Per-site pestPG (thetaStat print) has different columns — not used here
# Excludes DH (n=1) from collation
# =============================================================================
echo "Collating theta summary..."

Rscript --vanilla - "$THETA_DIR" << 'REOF'
args      <- commandArgs(trailingOnly = TRUE)
theta_dir <- args[1]

# Use windowed pestPG — has tW, tP, Tajima, nSites columns
# comment.char="" required — header line starts with #
# Exclude DH — n=1, population theta not meaningful
files <- list.files(theta_dir,
                    pattern    = "\\.thetas\\.windowed\\.pestPG$",
                    full.names = TRUE)
files <- files[!grepl("pop_map_DH", files)]

cat("Found", length(files), "windowed pestPG files\n")

results <- lapply(files, function(f) {
    pop <- sub("\\.thetas\\.windowed\\.pestPG$", "", basename(f))

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
        cat("WARNING: Empty or unreadable file:", f, "\n")
        return(NULL)
    }

    n_sites  <- sum(dat$nSites, na.rm = TRUE)
    theta_pi <- sum(dat$tP,     na.rm = TRUE)
    theta_w  <- sum(dat$tW,     na.rm = TRUE)
    tajima_D <- weighted.mean(dat$Tajima, dat$nSites, na.rm = TRUE)

    data.frame(
        population               = pop,
        n_sites                  = n_sites,
        Watterson_theta          = theta_w,
        theta_pi                 = theta_pi,
        tajima_D                 = tajima_D,
        Watterson_theta_per_site = theta_w  / n_sites,
        theta_pi_per_site        = theta_pi / n_sites,
        stringsAsFactors         = FALSE
    )
})

results <- results[!sapply(results, is.null)]

if (length(results) == 0) {
    cat("ERROR: No valid windowed pestPG files found\n")
    quit(status = 1)
}

theta_summary <- do.call(rbind, results)
theta_summary <- theta_summary[order(theta_summary$population), ]

out_file <- file.path(theta_dir, "theta_summary_corrected.tsv")
write.table(theta_summary, out_file,
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nTheta summary written to:", out_file, "\n\n")
print(theta_summary[, c("population", "n_sites",
                         "theta_pi_per_site",
                         "Watterson_theta_per_site",
                         "tajima_D")])
REOF

echo ""
echo "=========================================="
echo "Population theta rerun complete"
echo "Results in: $THETA_DIR"
echo "Summary: $THETA_DIR/theta_summary_corrected.tsv"
echo "=========================================="