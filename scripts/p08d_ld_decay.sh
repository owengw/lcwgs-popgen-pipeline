#!/bin/bash
#SBATCH --job-name=ld_decay
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/p08d_ld_decay_%j.log

# =============================================================================
# LD decay analysis using ngsLD
# Runs AFTER p08c_hwe.sh — reuses beagle GL files produced there
#
# Key requirements for this cluster:
#   - Use /users/bi4og/ngsLD/ngsLD (v1.2.1) NOT conda ngsLD (v1.2.0)
#   - Positions file must be gzipped and tab-separated (chr\tpos), NOT chr:pos
#   - LD_LIBRARY_PATH must include conda lib for libgsl.so.27 symlink
#   - libgsl.so.27 symlink must exist in conda lib dir:
#     ln -s .../libgsl.so.25 .../libgsl.so.27
#
# Populations: CAI, NS, SB, SI
# Distance: 100kb maximum
#
# POP_SINGLE variable supported:
#   sbatch --export=OUT_DIR=physalia/angsd_results,POP_SINGLE=SI p08d_ld_decay.sh
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

# Required for /users/bi4og/ngsLD/ngsLD binary — provides libgsl.so.27
export LD_LIBRARY_PATH=/mnt/parscratch/users/bi4og/conda_envs/owenspopgen/lib:${LD_LIBRARY_PATH:-}

module load R/4.4.1-foss-2022b 2>/dev/null || true

set -euo pipefail
set -x

# =============================================================================
# Variables
# =============================================================================
NGSLD="/users/bi4og/ngsLD/ngsLD"
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
REFERENCE="genome/Plong_genome_flye.fasta"
SITES_FILE="$OUT_DIR/random_10M_sites.list"
RF_FILE="$OUT_DIR/random_10M_contigs.txt"
POP_LIST_DIR="$OUT_DIR/pop_map_lists"
HWE_DIR="$OUT_DIR/hwe"
LD_DIR="$OUT_DIR/ld_decay"
THREADS=8

EXCLUDE="o_merged\|T_merged\|SI43_merged\|SI45_merged\|SI83_merged"

mkdir -p "$LD_DIR"
mkdir -p "$LD_DIR/clean_lists"
mkdir -p logs

# Validate ngsLD binary
if [[ ! -f "$NGSLD" ]]; then
    echo "ERROR: ngsLD binary not found: $NGSLD"
    exit 1
fi

# Ensure libgsl.so.27 symlink exists
GSL_SYMLINK="/mnt/parscratch/users/bi4og/conda_envs/owenspopgen/lib/libgsl.so.27"
GSL_SOURCE="/mnt/parscratch/users/bi4og/conda_envs/owenspopgen/lib/libgsl.so.25"
if [[ ! -f "$GSL_SYMLINK" ]]; then
    echo "Creating libgsl.so.27 symlink..."
    ln -s "$GSL_SOURCE" "$GSL_SYMLINK"
fi

POPULATIONS=(CAI NS SB SI)

if [[ -n "${POP_SINGLE:-}" ]]; then
    echo "POP_SINGLE=$POP_SINGLE — running single population only"
    POPULATIONS=("$POP_SINGLE")
fi

declare -A MIN_IND
MIN_IND[CAI]=7
MIN_IND[NS]=16
MIN_IND[SB]=5
MIN_IND[SI]=27

for POP in "${POPULATIONS[@]}"; do
    echo "=========================================="
    echo "Processing LD decay: $POP"
    echo "=========================================="

    FULL_LIST="$POP_LIST_DIR/pop_map_${POP}.list"
    CLEAN_LIST="$LD_DIR/clean_lists/${POP}_clean.list"
    BEAGLE_FILE="$HWE_DIR/pop_map_${POP}.beagle.gz"
    POS_FILE="$LD_DIR/${POP}_positions.gz"
    LD_OUT="$LD_DIR/${POP}_ld.tsv"
    DECAY_OUT="$LD_DIR/${POP}_ld_decay.tsv"

    if [[ ! -f "$FULL_LIST" ]]; then
        echo "WARNING: BAM list not found: $FULL_LIST — skipping"
        continue
    fi

    grep -v "$EXCLUDE" "$FULL_LIST" > "$CLEAN_LIST"
    N_CLEAN=$(wc -l < "$CLEAN_LIST")
    echo "  Clean samples: $N_CLEAN"

    # ------------------------------------------------------------------
    # Step 1 — Check beagle file from HWE run
    # ------------------------------------------------------------------
    if [[ ! -f "$BEAGLE_FILE" ]]; then
        echo "ERROR: Beagle file not found: $BEAGLE_FILE"
        echo "Run p08c_hwe.sh first"
        continue
    fi

    # Validate beagle file is not truncated
    N_SITES=$(zcat "$BEAGLE_FILE" | tail -n +2 | wc -l)
    if [[ $N_SITES -eq 0 ]]; then
        echo "ERROR: Beagle file is empty or truncated: $BEAGLE_FILE"
        continue
    fi
    echo "  Sites in beagle: $N_SITES"

    # ------------------------------------------------------------------
    # Step 2 — Generate positions file
    # Format: gzipped tab-separated chr\tpos (required by ngsLD v1.2.1)
    # Derived from beagle marker column: contig_XXX_POS -> contig_XXX\tPOS
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
    # --geno:        beagle GL file (gzipped)
    # --pos:         gzipped tab-separated positions (chr\tpos)
    # --probs:       input is genotype probabilities
    # --n_ind:       number of individuals
    # --n_sites:     number of sites
    # --max_kb_dist: maximum pairwise distance
    # --min_maf:     exclude rare variants
    # --rnd_sample:  fraction of pairs to compute (0.1 = 10%)
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
            tail -10 "$LD_DIR/${POP}_ngsld.log"
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

# ngsLD output columns: pos1 pos2 dist r^2 D D' r2_locus
if (grepl("\\.gz$", ld_file)) {
    dat <- read.table(gzfile(ld_file), header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
} else {
    dat <- read.table(ld_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE)
}

cat("Pairs loaded:", nrow(dat), "\n")
cat("Distance range:", min(dat$dist), "to", max(dat$dist), "bp\n")

# Bin by distance in 1kb windows up to 100kb
dat$dist_kb  <- dat$dist / 1000
breaks       <- seq(0, 100, by = 1)
labels       <- seq(0.5, 99.5, by = 1)
dat$dist_bin <- cut(dat$dist_kb, breaks = breaks,
                    labels = labels, include.lowest = TRUE)

decay <- do.call(rbind, lapply(levels(dat$dist_bin), function(b) {
    sub_dat <- dat[!is.na(dat$dist_bin) & dat$dist_bin == b, ]
    if (nrow(sub_dat) == 0) return(NULL)
    data.frame(
        population = pop,
        dist_kb    = as.numeric(b),
        mean_r2    = mean(sub_dat$r2,   na.rm = TRUE),
        median_r2  = median(sub_dat$r2, na.rm = TRUE),
        n_pairs    = nrow(sub_dat),
        stringsAsFactors = FALSE
    )
}))

decay <- decay[order(decay$dist_kb), ]

write.table(decay, out_file,
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("LD decay written to:", out_file, "\n")

# Report key distances
get_r2 <- function(target_kb) {
    idx <- which.min(abs(decay$dist_kb - target_kb))
    round(decay$mean_r2[idx], 4)
}
cat("Mean r2 at ~1kb: ",   get_r2(1),   "\n")
cat("Mean r2 at ~10kb:",   get_r2(10),  "\n")
cat("Mean r2 at ~50kb:",   get_r2(50),  "\n")
cat("Mean r2 at ~100kb:",  get_r2(99),  "\n")
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

files <- list.files(ld_dir, pattern = "_ld_decay\\.tsv$",
                    full.names = TRUE)

if (length(files) == 0) {
    cat("No decay files found\n")
    quit(status = 0)
}

combined <- do.call(rbind, lapply(files, function(f) {
    read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}))
combined <- combined[order(combined$population, combined$dist_kb), ]

out_file <- file.path(ld_dir, "ld_decay_all_populations.tsv")
write.table(combined, out_file,
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Combined LD decay written to:", out_file, "\n")

cat("\n=== LD at key distances (mean r2) ===\n")
for (pop in unique(combined$population)) {
    d <- combined[combined$population == pop, ]
    get_r2 <- function(kb) {
        round(d$mean_r2[which.min(abs(d$dist_kb - kb))], 4)
    }
    cat(pop, "— 1kb:", get_r2(1),
        "| 10kb:", get_r2(10),
        "| 50kb:", get_r2(50),
        "| 100kb:", get_r2(99), "\n")
}
REOF

echo ""
echo "=========================================="
echo "LD decay analysis complete"
echo "Results in: $LD_DIR"
echo "=========================================="
