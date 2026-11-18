#!/bin/bash
#SBATCH --job-name=hwe_ld
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/p08c_hwe_%j.log

# =============================================================================
# p08c_hwe.sh — RagTag rewrite.
# Generates per-population beagle GL files (random-sites subsample — HWE
# testing does not need every site, and this keeps the run tractable; LD
# decay in p08d reuses the same beagle), then tests HWE from GLs using the
# HardyWeinberg R package.
#
# CHANGES FROM ORIGINAL:
#  - REFERENCE is the real RagTag plongirostris.fasta, not Flye.
#  - No EXCLUDE regex. pop_map_lists/pop_map_{CAI,NS,SB,SI}.list already
#    reflect the final 180-sample set (SI35/replicates/dead libraries
#    already excluded, SI1T correctly kept).
#  - random_10M_sites.list/contigs regenerated fresh from chromosomes.txt
#    (the 13 placed RagTag chromosomes) if not already present — the old
#    file, if any, would be Flye-contig coordinates and unusable here.
#  - MIN_IND recomputed from ACTUAL current population sizes (CAI 22, NS 50,
#    SB 16, SI 85), not the stale repo defaults.
#
# POP_SINGLE supported: sbatch --export=OUT_DIR=...,POP_SINGLE=SI p08c_hwe.sh
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
module load R/4.4.1-foss-2022b 2>/dev/null || true
export R_LIBS_USER=/mnt/parscratch/users/bi4og/R_libs
set -eu
set -x

# =============================================================================
# Variables
# =============================================================================
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
REFERENCE="/mnt/parscratch/users/bi4og/ncbi_submission/plongirostris.fasta"
CHROM_FILE="$OUT_DIR/chromosomes.txt"
SITES_FILE="$OUT_DIR/random_10M_sites.list"
RF_FILE="$OUT_DIR/random_10M_contigs.txt"
POP_LIST_DIR="$OUT_DIR/pop_map_lists"
HWE_DIR="$OUT_DIR/hwe"
THREADS=${SLURM_CPUS_PER_TASK:-6}

mkdir -p "$HWE_DIR/clean_lists" logs

POPULATIONS=(CAI NS SB SI)
if [[ -n "${POP_SINGLE:-}" ]]; then
    echo "POP_SINGLE=$POP_SINGLE — running single population only"
    POPULATIONS=("$POP_SINGLE")
fi

# minInd = ~30% of each population's ACTUAL current size
declare -A MIN_IND
MIN_IND[CAI]=6    # 30% of 22
MIN_IND[NS]=12    # 30% of 50
MIN_IND[SB]=5     # 30% of 16
MIN_IND[SI]=25    # 30% of 85

# =============================================================================
# Regenerate random_10M_sites on RagTag if not already present
# (any pre-existing copy would be Flye-coordinate and unusable)
# =============================================================================
if [[ ! -f "$SITES_FILE" ]]; then
    echo "Generating random 10M sites from RagTag chromosomes >100kb..."
    [[ -f "$CHROM_FILE" ]] || { echo "ERROR: $CHROM_FILE not found"; exit 1; }
    [[ -f "$REFERENCE.fai" ]] || { echo "ERROR: $REFERENCE.fai not found"; exit 1; }
    awk 'NR==FNR{keep[$1]=1; next} ($1 in keep) && $2>100000' \
        "$CHROM_FILE" "$REFERENCE.fai" | \
        awk 'BEGIN{srand(42)}{for(i=1;i<=$2;i+=150) if(rand()<0.067) print $1"\t"i"\tA\tT"}' \
        | head -10000000 > "$SITES_FILE"
    angsd sites index "$SITES_FILE"
    cut -f1 "$SITES_FILE" | sort -u > "$RF_FILE"
    echo "  $(wc -l < "$SITES_FILE") sites written"
else
    echo "Using existing random sites file: $SITES_FILE ($(wc -l < "$SITES_FILE") sites)"
fi

# =============================================================================
# Install HardyWeinberg R package if not available (persistent user library)
# =============================================================================
R_USER_LIB="/mnt/parscratch/users/bi4og/R_libs"
mkdir -p "$R_USER_LIB"

Rscript --vanilla -e "
.libPaths(c('$R_USER_LIB', .libPaths()))
if (!requireNamespace('HardyWeinberg', quietly=TRUE)) {
    install.packages('HardyWeinberg', repos='https://cloud.r-project.org', lib='$R_USER_LIB')
    cat('HardyWeinberg installed\n')
} else {
    cat('HardyWeinberg already available\n')
}
"

for POP in "${POPULATIONS[@]}"; do
    echo "=========================================="
    echo "Processing HWE: $POP"
    echo "=========================================="

    FULL_LIST="$POP_LIST_DIR/pop_map_${POP}.list"
    OUT_PREFIX="$HWE_DIR/pop_map_${POP}"
    BEAGLE_FILE="${OUT_PREFIX}.beagle.gz"

    if [[ ! -f "$FULL_LIST" ]]; then
        echo "WARNING: BAM list not found: $FULL_LIST — skipping"
        continue
    fi

    N_CLEAN=$(wc -l < "$FULL_LIST")
    MIND=${MIN_IND[$POP]}
    echo "  Samples (current 180-set): $N_CLEAN, minInd: $MIND"

    # ------------------------------------------------------------------
    # Step 1 — Generate beagle GL file (shared with p08d LD decay)
    # ------------------------------------------------------------------
    if [[ -f "$BEAGLE_FILE" ]]; then
        N_SITES=$(zcat "$BEAGLE_FILE" | tail -n +2 | wc -l)
        FIRST_MARKER=$(zcat "$BEAGLE_FILE" | sed -n '2p' | cut -f1)
        if [[ "$FIRST_MARKER" != chr* ]]; then
            echo "  WARNING: existing beagle looks stale (marker: $FIRST_MARKER) — regenerating"
            rm -f "$BEAGLE_FILE"
        else
            echo "  Beagle already exists ($N_SITES sites, RagTag-confirmed), skipping ANGSD"
        fi
    fi

    if [[ ! -f "$BEAGLE_FILE" ]]; then
        angsd \
            -bam "$FULL_LIST" \
            -ref "$REFERENCE" \
            -GL 2 \
            -doGlf 2 \
            -doMajorMinor 1 \
            -doMaf 1 \
            -SNP_pval 1e-6 \
            -minMapQ 20 \
            -minQ 20 \
            -baq 1 \
            -C 50 \
            -uniqueOnly 1 \
            -minInd $MIND \
            -setMinDepthInd 1 \
            -remove_bads 1 \
            -trim 0 \
            -P $THREADS \
            -sites "$SITES_FILE" \
            -rf "$RF_FILE" \
            -out "$OUT_PREFIX" \
            2> "${OUT_PREFIX}.angsd.log"

        if [[ ! -f "$BEAGLE_FILE" ]]; then
            echo "ERROR: Beagle file not produced — check ${OUT_PREFIX}.angsd.log"
            tail -30 "${OUT_PREFIX}.angsd.log"
            continue
        fi
        N_SITES=$(zcat "$BEAGLE_FILE" | tail -n +2 | wc -l)
        echo "  Beagle sites: $N_SITES"
    fi

    N_CLEAN_FOR_R=$N_CLEAN

    # ------------------------------------------------------------------
    # Step 2 — HWE test from genotype likelihoods in R
    # ------------------------------------------------------------------
    echo "  Running HWE test in R..."

    R_LIBS_USER="$R_USER_LIB" Rscript --vanilla - "$BEAGLE_FILE" "$N_CLEAN_FOR_R" "$POP" "$HWE_DIR" << 'REOF'
args        <- commandArgs(trailingOnly = TRUE)
beagle_file <- args[1]
n_ind       <- as.integer(args[2])
pop         <- args[3]
hwe_dir     <- args[4]

user_lib <- Sys.getenv("R_LIBS_USER")
if (nchar(user_lib) > 0 && dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
if (!requireNamespace("HardyWeinberg", quietly = TRUE)) stop("HardyWeinberg not available")
library(HardyWeinberg)

cat("Loading beagle file:", beagle_file, "\n")
cat("Individuals (expected):", n_ind, "\n")

dat <- read.table(gzfile(beagle_file), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Sites loaded:", nrow(dat), "\n")

gl_cols  <- 4:ncol(dat)
n_gl_col <- length(gl_cols)
if (n_gl_col != n_ind * 3) {
    cat("WARNING: Expected", n_ind * 3, "GL columns, got", n_gl_col, "- adjusting n_ind\n")
    n_ind <- n_gl_col %/% 3
}
cat("Processing", nrow(dat), "sites across", n_ind, "individuals\n")

chunk_size <- 10000
n_sites    <- nrow(dat)
n_chunks   <- ceiling(n_sites / chunk_size)
cat("Processing in", n_chunks, "chunks of", chunk_size, "sites\n")

results <- vector("list", n_chunks)
for (chunk in seq_len(n_chunks)) {
    idx_start <- (chunk - 1) * chunk_size + 1
    idx_end   <- min(chunk * chunk_size, n_sites)
    chunk_dat <- dat[idx_start:idx_end, ]

    chunk_results <- lapply(seq_len(nrow(chunk_dat)), function(i) {
        gls <- as.numeric(chunk_dat[i, gl_cols])
        gl_mat <- matrix(gls, nrow = n_ind, ncol = 3, byrow = TRUE)
        row_sums <- rowSums(gl_mat); row_sums[row_sums == 0] <- 1
        gl_norm  <- gl_mat / row_sums

        n_AA <- sum(gl_norm[, 1], na.rm = TRUE)
        n_Aa <- sum(gl_norm[, 2], na.rm = TRUE)
        n_aa <- sum(gl_norm[, 3], na.rm = TRUE)
        counts <- round(c(AA = n_AA, AB = n_Aa, BB = n_aa))
        if (counts["AB"] == 0 && (counts["AA"] == 0 || counts["BB"] == 0)) return(NULL)

        p_val <- tryCatch(HWExact(counts, alternative = "two.sided", verbose = FALSE)$pval,
                          error = function(e) NA_real_)
        p_A  <- (2 * counts["AA"] + counts["AB"]) / (2 * sum(counts))
        p_a  <- 1 - p_A
        H_exp <- 2 * p_A * p_a
        H_obs <- counts["AB"] / sum(counts)
        F_site <- ifelse(H_exp > 0, 1 - H_obs / H_exp, NA_real_)

        data.frame(marker = chunk_dat$marker[i], n_AA = counts["AA"], n_Aa = counts["AB"],
                   n_aa = counts["BB"], H_obs = round(H_obs, 4), H_exp = round(H_exp, 4),
                   F_site = round(F_site, 4), pHWE = p_val,
                   stringsAsFactors = FALSE, row.names = NULL)
    })
    results[[chunk]] <- do.call(rbind, chunk_results[!sapply(chunk_results, is.null)])
    if (chunk %% 10 == 0) cat("  Processed chunk", chunk, "/", n_chunks, "\n")
}

hwe_results <- do.call(rbind, results[!sapply(results, is.null)])
cat("Sites with HWE results:", nrow(hwe_results), "\n")

hwe_file <- file.path(hwe_dir, paste0("pop_map_", pop, "_hwe_results.tsv.gz"))
gz_con   <- gzfile(hwe_file, "w")
write.table(hwe_results, gz_con, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz_con)
cat("Per-site HWE results written to:", hwe_file, "\n")

n_tested   <- nrow(hwe_results)
bonf_thresh <- 0.05 / n_tested
n_sig_005  <- sum(hwe_results$pHWE < 0.05,  na.rm = TRUE)
n_sig_001  <- sum(hwe_results$pHWE < 0.001, na.rm = TRUE)
n_sig_bonf <- sum(hwe_results$pHWE < bonf_thresh, na.rm = TRUE)

cat("\n=== HWE Summary:", pop, "===\n")
cat("Sites tested:          ", n_tested, "\n")
cat("Significant p<0.05:    ", n_sig_005, "(", round(100*n_sig_005/n_tested,2), "%)\n")
cat("Significant p<0.001:   ", n_sig_001, "\n")
cat("Significant Bonferroni:", n_sig_bonf, "(threshold:", round(bonf_thresh,8), ")\n")
cat("Mean per-site F:       ", round(mean(hwe_results$F_site, na.rm=TRUE), 4), "\n")
cat("Median per-site F:     ", round(median(hwe_results$F_site, na.rm=TRUE), 4), "\n")

summary_df <- data.frame(population = pop, n_sites_tested = n_tested,
    n_sig_p005 = n_sig_005, n_sig_p001 = n_sig_001, n_sig_bonferroni = n_sig_bonf,
    bonferroni_threshold = bonf_thresh, pct_sig_p005 = round(100*n_sig_005/n_tested, 2),
    mean_F_persite = round(mean(hwe_results$F_site, na.rm=TRUE), 6),
    median_F_persite = round(median(hwe_results$F_site, na.rm=TRUE), 6),
    stringsAsFactors = FALSE)

write.table(summary_df, file.path(hwe_dir, paste0(pop, "_hwe_summary.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Summary written\n")
REOF

    echo "  Done: $POP"
done

# =============================================================================
# Collate all population HWE summaries
# =============================================================================
echo ""
echo "Collating HWE summaries..."

Rscript --vanilla - "$HWE_DIR" << 'REOF'
args    <- commandArgs(trailingOnly = TRUE)
hwe_dir <- args[1]
files <- list.files(hwe_dir, pattern = "_hwe_summary\\.tsv$", full.names = TRUE)
if (length(files) == 0) { cat("No HWE summary files found\n"); quit(status = 0) }
combined <- do.call(rbind, lapply(files, function(f) read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)))
combined <- combined[order(combined$population), ]
out_file <- file.path(hwe_dir, "hwe_summary_all_populations.tsv")
write.table(combined, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n=== HWE Summary across populations ===\n"); print(combined)
cat("\nWritten to:", out_file, "\n")
REOF

echo ""
echo "=========================================="
echo "HWE analysis complete"
echo "Results in: $HWE_DIR"
echo "Beagle files also available for LD decay: $HWE_DIR/pop_map_*.beagle.gz"
echo "=========================================="
