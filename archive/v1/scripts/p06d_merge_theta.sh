#!/bin/bash
#SBATCH --job-name=p06d_merge_theta
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/p06d_merge_theta_%A.log

# Load environment
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -euo pipefail
set -x

# Variables with defaults
BAM_LIST=${BAM_LIST:-}
REFERENCE=${REFERENCE:-$PWD/genome/Plong_genome_flye.fasta}
OUT_DIR=${OUT_DIR:-$PWD/physalia/angsd_results}
ANCESTRAL=${ANCESTRAL:-}
MIN_MAPQ=${MIN_MAPQ:-20}
MIN_Q=${MIN_Q:-20}
FOLD_SFS=${FOLD_SFS:-1}
WINDOW_SIZE=${WINDOW_SIZE:-50000}
WINDOW_STEP=${WINDOW_STEP:-10000}

PREFIX=$(basename "$BAM_LIST" .list)
OUT_PREFIX="$OUT_DIR/$PREFIX"

echo "Starting theta calculation for $PREFIX"

# Ensure binaries are available
REALSFS_BIN=${REALSFS_BIN:-realSFS}
THETASTAT_BIN=${THETASTAT_BIN:-thetaStat}

# -------------------------------
# Validate inputs
# -------------------------------
if [[ ! -f "$OUT_PREFIX.saf.idx" ]]; then
    echo "ERROR: SAF not found: $OUT_PREFIX.saf.idx"
    echo "Run p06c_angsd_subset.sh first for $PREFIX"
    exit 1
fi

# -------------------------------
# Estimate SFS (if not already done)
# -------------------------------
if [[ ! -f "$OUT_PREFIX.sfs" ]]; then
    echo "Estimating SFS..."
    if ! $REALSFS_BIN "$OUT_PREFIX.saf.idx" \
        -fold $FOLD_SFS \
        -P 4 > "$OUT_PREFIX.sfs" 2> "$OUT_PREFIX.sfs.log"; then
        echo "ERROR: SFS estimation failed"
        cat "$OUT_PREFIX.sfs.log"
        exit 1
    fi
    echo "? SFS estimated"
else
    echo "SFS already exists, using existing file"
fi

# -------------------------------
# Calculate thetas using realSFS saf2theta
# -------------------------------
echo "Running realSFS saf2theta..."

# Build ancestral option
ANC_OPT=""
if [[ -n "$ANCESTRAL" && -f "$ANCESTRAL" ]]; then
    ANC_OPT="-anc $ANCESTRAL"
    echo "Using ancestral genome for polarization"
else
    echo "No ancestral genome, using folded spectrum"
fi

# Run saf2theta
if [[ -n "$ANC_OPT" ]]; then
    if ! $REALSFS_BIN saf2theta "$OUT_PREFIX.saf.idx" \
        -sfs "$OUT_PREFIX.sfs" \
        -outname "$OUT_PREFIX" \
        -P 4 \
        $ANC_OPT > "$OUT_PREFIX.saf2theta.log" 2>&1; then
        echo "ERROR: Failed to calculate thetas with saf2theta"
        cat "$OUT_PREFIX.saf2theta.log"
        exit 1
    fi
else
    if ! $REALSFS_BIN saf2theta "$OUT_PREFIX.saf.idx" \
        -sfs "$OUT_PREFIX.sfs" \
        -outname "$OUT_PREFIX" \
        -P 4 \
        -fold $FOLD_SFS > "$OUT_PREFIX.saf2theta.log" 2>&1; then
        echo "ERROR: Failed to calculate thetas with saf2theta"
        cat "$OUT_PREFIX.saf2theta.log"
        exit 1
    fi
fi

# -------------------------------
# Verify theta index created
# -------------------------------
if [[ ! -f "$OUT_PREFIX.thetas.idx" ]]; then
    echo "ERROR: Theta index file not created: $OUT_PREFIX.thetas.idx"
    exit 1
fi

echo "? Thetas calculated"

# -------------------------------
# Theta statistics
# -------------------------------
echo "Calculating theta statistics..."

# Genome-wide statistics
if ! $THETASTAT_BIN do_stat "$OUT_PREFIX.thetas.idx" > "$OUT_PREFIX.thetas.stats.log" 2>&1; then
    echo "ERROR: Failed to calculate genome-wide theta statistics"
    cat "$OUT_PREFIX.thetas.stats.log"
    exit 1
fi

# Windowed statistics
if ! $THETASTAT_BIN do_stat "$OUT_PREFIX.thetas.idx" \
    -win $WINDOW_SIZE -step $WINDOW_STEP \
    -outnames "$OUT_PREFIX.thetas.windowed" > "$OUT_PREFIX.thetas.windowed.log" 2>&1; then
    echo "ERROR: Failed to calculate windowed theta statistics"
    cat "$OUT_PREFIX.thetas.windowed.log"
    exit 1
fi

# thetaStat appends .pestPG to the output name
if [[ -f "$OUT_PREFIX.thetas.windowed.pestPG" ]]; then
    mv "$OUT_PREFIX.thetas.windowed.pestPG" "$OUT_PREFIX.thetas.windowed.tsv"
    echo "? Windowed theta statistics calculated"
else
    echo "ERROR: Expected output file not found: $OUT_PREFIX.thetas.windowed.pestPG"
    exit 1
fi

echo ""
echo "=================================================================="
echo "Theta calculation complete for $PREFIX"
echo "=================================================================="
echo "Output files:"
echo "  - $OUT_PREFIX.thetas.idx"
echo "  - $OUT_PREFIX.thetas.idx.pestPG (genome-wide)"
echo "  - $OUT_PREFIX.thetas.windowed.tsv (windowed)"
echo "=================================================================="