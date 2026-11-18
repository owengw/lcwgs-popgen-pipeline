#!/bin/bash
#SBATCH --job-name=hwe_beagle_chr
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --output=logs/hwe_beagle_chr_%A_%a.log

# =============================================================================
# p08c_hwe_bychrom.sh — per-chromosome HWE beagle generation, avoiding the
# no-checkpoint problem: ANGSD's serial genome-wide scan means any timeout
# forces a full restart from chr1. Splitting by chromosome (-r CHR, matching
# the -sites/-rf restriction already used for random_10M_sites) means each
# chromosome only needs redoing once if it individually times out.
#
# Modes:
#   --array-task   (ARRAY TASK) generate one chromosome's beagle chunk
#   --concat       concatenate all 13 chunks -> final population beagle,
#                  then run the same R HWE test as p08c_hwe.sh
#
# USAGE:
#   NCHR=$(wc -l < chromosomes.txt)
#   sbatch --array=1-$NCHR%13 --export=ALL,POP=SI OUT_DIR=... p08c_hwe_bychrom.sh --array-task
#   sbatch --wrap="bash p08c_hwe_bychrom.sh --concat POP=SI OUT_DIR=..."
# =============================================================================

if [[ "${1:-}" == "--array-task" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    : "${POP:?POP required (e.g. SI)}"
    source ~/.bash_profile
    conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
    set -eu
    THREADS=${SLURM_CPUS_PER_TASK:-6}
    OUT_DIR=${OUT_DIR:-physalia/angsd_results}
    REFERENCE=${REFERENCE:-/mnt/parscratch/users/bi4og/ncbi_submission/plongirostris.fasta}
    CHROM_FILE=${CHROM_FILE:-$OUT_DIR/chromosomes.txt}
    SITES_FILE=${SITES_FILE:-$OUT_DIR/random_10M_sites.list}
    FULL_LIST=${FULL_LIST:-$OUT_DIR/pop_map_lists/pop_map_${POP}.list}
    MIND=${MIND:?MIND required (e.g. 26 for SI)}

    CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_FILE")
    [[ -n "$CHR" ]] || { echo "ERROR: no chromosome on line $SLURM_ARRAY_TASK_ID"; exit 1; }

    mkdir -p "$OUT_DIR/hwe/beagle_parts"
    OUT_PREFIX="$OUT_DIR/hwe/beagle_parts/${POP}.${CHR}"

    if [[ -s "${OUT_PREFIX}.beagle.gz" ]]; then
        # verify it's not a truncated leftover before trusting the skip
        if gzip -t "${OUT_PREFIX}.beagle.gz" 2>/dev/null; then
            echo "Already done: ${OUT_PREFIX}.beagle.gz"; exit 0
        else
            echo "Existing beagle is truncated, regenerating: ${OUT_PREFIX}.beagle.gz"
            rm -f "${OUT_PREFIX}.beagle.gz"
        fi
    fi

    echo "=================================================================="
    echo "HWE beagle  pop=$POP  chr=$CHR  minInd=$MIND  threads=$THREADS"
    echo "=================================================================="

    angsd -bam "$FULL_LIST" -ref "$REFERENCE" \
        -r "$CHR" \
        -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 \
        -SNP_pval 1e-6 \
        -minMapQ 20 -minQ 20 -baq 1 -C 50 -uniqueOnly 1 \
        -minInd "$MIND" -setMinDepthInd 1 -remove_bads 1 -trim 0 \
        -P "$THREADS" \
        -sites "$SITES_FILE" \
        -out "$OUT_PREFIX" \
        2> "${OUT_PREFIX}.angsd.log"

    gzip -t "${OUT_PREFIX}.beagle.gz" || { echo "ERROR: beagle for $CHR is corrupt after generation"; exit 1; }
    echo "Done: ${OUT_PREFIX}.beagle.gz"
    exit 0
fi

if [[ "${1:-}" == "--concat" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    : "${POP:?POP required}"
    source ~/.bash_profile
    conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
    module load R/4.4.1-foss-2022b 2>/dev/null || true
    export R_LIBS_USER=""
    set -eu
    OUT_DIR=${OUT_DIR:-physalia/angsd_results}
    CHROM_FILE=${CHROM_FILE:-$OUT_DIR/chromosomes.txt}
    HWE_DIR="$OUT_DIR/hwe"
    N_CLEAN=${N_CLEAN:?N_CLEAN required (individual count for this population)}

    PARTS=()
    while read -r c; do
        f="$HWE_DIR/beagle_parts/${POP}.${c}.beagle.gz"
        [[ -s "$f" ]] || { echo "ERROR: missing chromosome part: $f"; exit 1; }
        gzip -t "$f" || { echo "ERROR: corrupt part: $f"; exit 1; }
        PARTS+=("$f")
    done < "$CHROM_FILE"
    echo "All ${#PARTS[@]} chromosome parts present and valid"

    BEAGLE_FILE="$HWE_DIR/pop_map_${POP}.beagle.gz"
    { zcat "${PARTS[0]}" | head -1
      for f in "${PARTS[@]}"; do zcat "$f" | tail -n +2; done
    } | gzip > "$BEAGLE_FILE"

    N_SITES=$(zcat "$BEAGLE_FILE" | tail -n +2 | wc -l)
    echo "Concatenated beagle: $N_SITES sites -> $BEAGLE_FILE"
    [[ "$N_SITES" -gt 0 ]] || { echo "ERROR: concatenated beagle is empty"; exit 1; }

    # ---- same HWE R test as p08c_hwe.sh ----
    R_LIBS_USER="$R_LIBS_USER" Rscript --vanilla - "$BEAGLE_FILE" "$N_CLEAN" "$POP" "$HWE_DIR" << 'REOF'
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
dat <- read.table(gzfile(beagle_file), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Sites loaded:", nrow(dat), "\n")

gl_cols  <- 4:ncol(dat)
n_gl_col <- length(gl_cols)
if (n_gl_col != n_ind * 3) {
    cat("WARNING: Expected", n_ind * 3, "GL columns, got", n_gl_col, "- adjusting n_ind\n")
    n_ind <- n_gl_col %/% 3
}

chunk_size <- 10000
n_sites    <- nrow(dat)
n_chunks   <- ceiling(n_sites / chunk_size)

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
}

hwe_results <- do.call(rbind, results[!sapply(results, is.null)])
cat("Sites with HWE results:", nrow(hwe_results), "\n")

hwe_file <- file.path(hwe_dir, paste0("pop_map_", pop, "_hwe_results.tsv.gz"))
gz_con   <- gzfile(hwe_file, "w")
write.table(hwe_results, gz_con, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz_con)

n_tested   <- nrow(hwe_results)
bonf_thresh <- 0.05 / n_tested
n_sig_005  <- sum(hwe_results$pHWE < 0.05,  na.rm = TRUE)
n_sig_001  <- sum(hwe_results$pHWE < 0.001, na.rm = TRUE)
n_sig_bonf <- sum(hwe_results$pHWE < bonf_thresh, na.rm = TRUE)

summary_df <- data.frame(population = pop, n_sites_tested = n_tested,
    n_sig_p005 = n_sig_005, n_sig_p001 = n_sig_001, n_sig_bonferroni = n_sig_bonf,
    bonferroni_threshold = bonf_thresh, pct_sig_p005 = round(100*n_sig_005/n_tested, 2),
    mean_F_persite = round(mean(hwe_results$F_site, na.rm=TRUE), 6),
    median_F_persite = round(median(hwe_results$F_site, na.rm=TRUE), 6),
    stringsAsFactors = FALSE)

write.table(summary_df, file.path(hwe_dir, paste0(pop, "_hwe_summary.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Summary written\n"); print(summary_df)
REOF

    echo "Done: $POP"
    exit 0
fi

echo "Usage: $0 --array-task|--concat POP=<pop> OUT_DIR=<dir> [MIND=<n>] [N_CLEAN=<n>]"
exit 1
