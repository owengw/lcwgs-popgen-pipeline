#!/bin/bash
#SBATCH --job-name=ld_decay
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/p08d_ld_decay_%j.log

# =============================================================================
# p08d_ld_decay.sh — RagTag rewrite.
# Reuses the per-population beagle GL files produced by p08c_hwe.sh
# ($OUT_DIR/hwe/pop_map_<POP>.beagle.gz) — no new ANGSD run.
#
# CHANGES FROM ORIGINAL:
#  - No EXCLUDE regex. pop_map_lists/pop_map_{CAI,NS,SB,SI}.list already
#    reflect the final 180-sample set.
#  - N_CLEAN read directly from the pop_map list (no separate clean-list
#    rebuild step — that filtering already happened at the sample-list
#    stage this session, not per-script).
#  - Stale-beagle guard: if the beagle's first marker isn't chr*, refuses
#    to proceed rather than silently computing LD decay on Flye coordinates.
#
# Populations: CAI, NS, SB, SI. Distance: 100kb maximum.
# POP_SINGLE supported: sbatch --export=OUT_DIR=...,POP_SINGLE=SI p08d_ld_decay.sh
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

export LD_LIBRARY_PATH=/mnt/parscratch/users/bi4og/conda_envs/owenspopgen/lib:${LD_LIBRARY_PATH:-}
module load R/4.4.1-foss-2022b 2>/dev/null || true
export R_LIBS_USER=/mnt/parscratch/users/bi4og/R_libs
set -eu
set -x

NGSLD="/users/bi4og/ngsLD/ngsLD"
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
POP_LIST_DIR="$OUT_DIR/pop_map_lists"
HWE_DIR="$OUT_DIR/hwe"
LD_DIR="$OUT_DIR/ld_decay"
THREADS=8

mkdir -p "$LD_DIR" logs

[[ -f "$NGSLD" ]] || { echo "ERROR: ngsLD binary not found: $NGSLD"; exit 1; }

GSL_SYMLINK="/mnt/parscratch/users/bi4og/conda_envs/owenspopgen/lib/libgsl.so.27"
GSL_SOURCE="/mnt/parscratch/users/bi4og/conda_envs/owenspopgen/lib/libgsl.so.25"
if [[ ! -f "$GSL_SYMLINK" && -f "$GSL_SOURCE" ]]; then
    echo "Creating libgsl.so.27 symlink..."
    ln -s "$GSL_SOURCE" "$GSL_SYMLINK"
fi

POPULATIONS=(CAI NS SB SI)
if [[ -n "${POP_SINGLE:-}" ]]; then
    echo "POP_SINGLE=$POP_SINGLE — running single population only"
    POPULATIONS=("$POP_SINGLE")
fi

for POP in "${POPULATIONS[@]}"; do
    echo "=========================================="
    echo "Processing LD decay: $POP"
    echo "=========================================="

    FULL_LIST="$POP_LIST_DIR/pop_map_${POP}.list"
    BEAGLE_FILE="$HWE_DIR/pop_map_${POP}.beagle.gz"
    POS_FILE="$LD_DIR/${POP}_positions.gz"
    LD_OUT="$LD_DIR/${POP}_ld.tsv"
    DECAY_OUT="$LD_DIR/${POP}_ld_decay.tsv"

    if [[ ! -f "$FULL_LIST" ]]; then
        echo "WARNING: BAM list not found: $FULL_LIST — skipping"
        continue
    fi
    N_CLEAN=$(wc -l < "$FULL_LIST")
    echo "  Samples (current 180-set): $N_CLEAN"

    # ------------------------------------------------------------------
    # Step 1 — Beagle must already exist (from p08c) AND be RagTag, not stale
    # ------------------------------------------------------------------
    if [[ ! -f "$BEAGLE_FILE" ]]; then
        echo "ERROR: Beagle file not found: $BEAGLE_FILE"
        echo "Run p08c_hwe.sh first"
        continue
    fi

    FIRST_MARKER=$(zcat "$BEAGLE_FILE" | sed -n '2p' | cut -f1)
    if [[ "$FIRST_MARKER" != chr* ]]; then
        echo "ERROR: Beagle looks stale (marker: $FIRST_MARKER, expected chr*_RagTag_pilon...)"
        echo "Rerun p08c_hwe.sh to regenerate it before continuing"
        continue
    fi

    N_SITES=$(zcat "$BEAGLE_FILE" | tail -n +2 | wc -l)
    if [[ $N_SITES -eq 0 ]]; then
        echo "ERROR: Beagle file is empty: $BEAGLE_FILE"
        continue
    fi
    echo "  Sites in beagle: $N_SITES (RagTag-confirmed: $FIRST_MARKER)"

    # ------------------------------------------------------------------
    # Step 2 — Generate positions file (gzipped tab-separated chr\tpos)
    # ------------------------------------------------------------------
    if [[ -f "$POS_FILE" ]]; then
        echo "  Positions file already exists, skipping"
    else
        echo "  Generating positions file..."
        zcat "$BEAGLE_FILE" | tail -n +2 | \
            awk '{
                n = split($1, a, "_")
                chr = ""
                for (i = 1; i < n; i++) chr = chr (i > 1 ? "_" : "") a[i]
                print chr "\t" a[n]
            }' | gzip > "$POS_FILE"

        N_POS=$(zcat "$POS_FILE" | wc -l)
        echo "  Positions: $N_POS (should match $N_SITES)"
        if [[ $N_POS -ne $N_SITES ]]; then
            echo "ERROR: Position count mismatch — $N_POS vs $N_SITES sites"
            exit 1
        fi
    fi

    # ------------------------------------------------------------------
    # Step 3 — Run ngsLD
    # ------------------------------------------------------------------
    if [[ -f "$LD_OUT" ]] || [[ -f "${LD_OUT}.gz" ]]; then
        echo "  ngsLD output already exists, skipping"
    else
        echo "  Running ngsLD..."
        "$NGSLD" \
            --geno "$BEAGLE_FILE" \
            --pos "$POS_FILE" \
            --probs \
            --n_ind $N_CLEAN \
            --n_sites $N_SITES \
            --max_kb_dist 100 \
            --min_maf 0.05 \
            --rnd_sample 0.1 \
            --n_threads $THREADS \
            --out "$LD_OUT" \
            2> "$LD_DIR/${POP}_ngsld.log"

        if [[ ! -f "$LD_OUT" ]]; then
            echo "ERROR: ngsLD failed — check $LD_DIR/${POP}_ngsld.log"
            tail -20 "$LD_DIR/${POP}_ngsld.log"
            continue
        fi

        N_PAIRS=$(wc -l < "$LD_OUT")
        echo "  LD pairs: $N_PAIRS"
        gzip "$LD_OUT"
        echo "  Output gzipped: ${LD_OUT}.gz"
    fi

    # ------------------------------------------------------------------
    # Step 4 — Bin LD by distance and compute decay curve
    # ------------------------------------------------------------------
    echo "  Computing LD decay curve..."
    LD_FILE="${LD_OUT}.gz"
    [[ ! -f "$LD_FILE" ]] && LD_FILE="$LD_OUT"

    Rscript --vanilla - "$LD_FILE" "$DECAY_OUT" "$POP" << 'REOF'
args     <- commandArgs(trailingOnly = TRUE)
ld_file  <- args[1]
out_file <- args[2]
pop      <- args[3]

cat("Loading LD data for:", pop, "\n")
con <- if (grepl("\\.gz$", ld_file)) gzfile(ld_file) else ld_file
dat <- read.table(con, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                  col.names = c("pos1","pos2","dist","r2","D","Dprime","r2_expected"))

cat("Pairs loaded:", nrow(dat), "\n")
cat("Distance range:", min(dat$dist, na.rm=TRUE), "to", max(dat$dist, na.rm=TRUE), "bp\n")

dat$dist_kb  <- dat$dist / 1000
breaks       <- seq(0, 100, by = 1)
labels       <- seq(0.5, 99.5, by = 1)
dat$dist_bin <- cut(dat$dist_kb, breaks = breaks, labels = labels, include.lowest = TRUE)

decay <- do.call(rbind, lapply(levels(dat$dist_bin), function(b) {
    sub_dat <- dat[!is.na(dat$dist_bin) & dat$dist_bin == b, ]
    if (nrow(sub_dat) == 0) return(NULL)
    data.frame(population = pop, dist_kb = as.numeric(b),
               mean_r2 = mean(sub_dat$r2, na.rm = TRUE),
               median_r2 = median(sub_dat$r2, na.rm = TRUE),
               n_pairs = nrow(sub_dat), stringsAsFactors = FALSE)
}))
decay <- decay[order(decay$dist_kb), ]
write.table(decay, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("LD decay written to:", out_file, "\n")

get_r2 <- function(target_kb) {
    idx <- which.min(abs(decay$dist_kb - target_kb))
    round(decay$mean_r2[idx], 4)
}
cat("Mean r2 at ~1kb: ",  get_r2(1),  "\n")
cat("Mean r2 at ~10kb:",  get_r2(10), "\n")
cat("Mean r2 at ~50kb:",  get_r2(50), "\n")
cat("Mean r2 at ~100kb:", get_r2(99), "\n")
REOF

    echo "  Done: $POP"
done

# =============================================================================
# Combine decay curves across populations
# =============================================================================
echo ""
echo "Combining LD decay curves..."

Rscript --vanilla - "$LD_DIR" << 'REOF'
args   <- commandArgs(trailingOnly = TRUE)
ld_dir <- args[1]
files <- list.files(ld_dir, pattern = "_ld_decay\\.tsv$", full.names = TRUE)
if (length(files) == 0) { cat("No decay files found\n"); quit(status = 0) }
combined <- do.call(rbind, lapply(files, function(f) read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)))
combined <- combined[order(combined$population, combined$dist_kb), ]
out_file <- file.path(ld_dir, "ld_decay_all_populations.tsv")
write.table(combined, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Combined LD decay written to:", out_file, "\n")

cat("\n=== LD at key distances (mean r2) ===\n")
for (pop in unique(combined$population)) {
    d <- combined[combined$population == pop, ]
    get_r2 <- function(kb) round(d$mean_r2[which.min(abs(d$dist_kb - kb))], 4)
    cat(pop, "- 1kb:", get_r2(1), "| 10kb:", get_r2(10), "| 50kb:", get_r2(50), "| 100kb:", get_r2(99), "\n")
}
REOF

echo ""
echo "=========================================="
echo "LD decay analysis complete"
echo "Results in: $LD_DIR"
echo "=========================================="
