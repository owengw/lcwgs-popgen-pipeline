#!/bin/bash
#SBATCH --job-name=p06c_angsd
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p06c_angsd_%A.log

# Load environment
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen || true

set -euo pipefail
set -x

# Define binaries
ANGSD_BIN=angsd

# Variables
BAM_LIST=${BAM_LIST:-}
REFERENCE=${REFERENCE:-}
OUT_DIR=${OUT_DIR:-}
ANCESTRAL=${ANCESTRAL:-}
THREADS=${THREADS:-8}
MIN_MAPQ=${MIN_MAPQ:-10}
MIN_Q=${MIN_Q:-13}
FOLD_SFS=${FOLD_SFS:-1}
MIN_IND_RATIO=${MIN_IND_RATIO:-0.20}
USE_RANDOM_SITES=${USE_RANDOM_SITES:-0}  # Optional random site subsampling

# Validate inputs
if [[ ! -f "$BAM_LIST" ]]; then
    echo "ERROR: BAM list not found: $BAM_LIST"
    exit 1
fi

if [[ ! -f "$REFERENCE.fai" ]]; then
    echo "ERROR: Reference FAI not found at $REFERENCE.fai"
    exit 1
fi

# Extract population/group name
POP_NAME=$(basename "$BAM_LIST" .list)
OUT_PREFIX="$OUT_DIR/$POP_NAME"

mkdir -p "$OUT_DIR/logs"

# Calculate minimum individuals
N_SAMPLES=$(wc -l < "$BAM_LIST")
MIN_IND=$(echo "$N_SAMPLES * $MIN_IND_RATIO" | bc | awk '{print int($1+0.5)}')

# Sanity checks
[[ $MIN_IND -lt 1 ]] && MIN_IND=1
[[ $MIN_IND -gt $N_SAMPLES ]] && MIN_IND=$N_SAMPLES

echo "=================================================================="
echo "ANGSD: $POP_NAME (SAF Generation)"
echo "=================================================================="
echo "Population: $POP_NAME"
echo "Samples: $N_SAMPLES"
echo "Min individuals: $MIN_IND"
echo "Threads: $THREADS"
echo "=================================================================="

# -------------------------------
# Ancestral genome setup
# -------------------------------
if [[ -n "$ANCESTRAL" && -f "$ANCESTRAL" ]]; then
    ANC_OPT="-anc $ANCESTRAL"
    echo "Using ancestral genome: $ANCESTRAL"
else
    ANC_OPT="-anc $REFERENCE"
    echo "Using reference as ancestral (folded): $REFERENCE"
fi

# -------------------------------
# Random site subsampling setup
# -------------------------------
SITES_OPT=""
RF_OPT=""

if [[ "$USE_RANDOM_SITES" == "1" ]]; then
    # Generate random sites file if it doesn't exist
    if [[ ! -f "$OUT_DIR/random_10M_sites.list" ]]; then
        echo "Generating memory-efficient random 10M site subsample..."
        echo "Strategy: Sampling from largest contigs only to reduce memory usage"
        
        # Step 1: Generate random sites from LARGEST contigs only (>100kb)
        # This concentrates sites in fewer contigs, reducing memory dramatically
        awk '$2 > 100000' "$REFERENCE.fai" | \
            awk 'BEGIN{srand(42)} {
                # Sample uniformly: every ~150bp, 6.7% chance
                for(i=1; i<=$2; i+=150) {
                    if(rand() < 0.067) print $1"\t"i"\tA\tT"
                }
            }' | head -10000000 > "$OUT_DIR/random_10M_sites.list"
        
        echo "Generated $(wc -l < $OUT_DIR/random_10M_sites.list) random sites"
        
        # Step 2: Index for ANGSD
        echo "Indexing sites file..."
        angsd sites index "$OUT_DIR/random_10M_sites.list"
        
        # Step 3: Create contig list for -rf flag (CRITICAL for memory reduction)
        echo "Creating contig restriction list..."
        cut -f1 "$OUT_DIR/random_10M_sites.list" | sort -u > "$OUT_DIR/random_10M_contigs.txt"
        
        echo "Sites span $(wc -l < $OUT_DIR/random_10M_contigs.txt) contigs"
        
    else
        echo "Using existing random sites file"
        
        # Ensure contig list exists
        if [[ ! -f "$OUT_DIR/random_10M_contigs.txt" ]]; then
            echo "Regenerating contig list from existing sites..."
            cut -f1 "$OUT_DIR/random_10M_sites.list" | sort -u > "$OUT_DIR/random_10M_contigs.txt"
        fi
    fi
    
    SITES_OPT="-sites $OUT_DIR/random_10M_sites.list"
    RF_OPT="-rf $OUT_DIR/random_10M_contigs.txt"
    
    echo "Using random 10M site subsample (genome-wide, unbiased)"
    echo "Memory optimization: -rf flag limits ANGSD to $(wc -l < $OUT_DIR/random_10M_contigs.txt) contigs"
fi

# -------------------------------
# ANGSD Run: SAF generation ONLY
# Genome-wide (or random subset) for unbiased theta/FST
# -------------------------------
echo ""
echo "Running ANGSD (SAF generation for theta/FST)..."

if ! $ANGSD_BIN -bam "$BAM_LIST" -ref "$REFERENCE" $ANC_OPT \
    -GL 2 -doSaf 1 -doCounts 1 \
    -minMapQ $MIN_MAPQ -minQ $MIN_Q -minInd $MIN_IND \
    -setMinDepthInd 1 -P $THREADS \
    -nQueueSize 50 -remove_bads 1 -trim 0 \
    $SITES_OPT \
    $RF_OPT \
    -out "$OUT_PREFIX" 2>&1 | tee "$OUT_PREFIX.angsd_saf.log"; then
    echo "ERROR: ANGSD SAF generation failed"
    tail -50 "$OUT_PREFIX.angsd_saf.log"
    exit 1
fi

# -------------------------------
# Verify SAF files
# -------------------------------
if [[ ! -f "$OUT_PREFIX.saf.idx" ]]; then
    echo "ERROR: SAF index not created"
    exit 1
fi

echo "? SAF generation completed successfully"

# Show summary
echo ""
echo "SAF summary:"
grep -E "sites|Total number" "$OUT_PREFIX.angsd_saf.log" | tail -10

echo ""
echo "=================================================================="
echo "ANGSD SAF Complete: $POP_NAME"
echo "=================================================================="
echo "Output files:"
echo "  SAF index:  $OUT_PREFIX.saf.idx"
echo "  SAF values: $OUT_PREFIX.saf.gz"
echo "  SAF pos:    $OUT_PREFIX.saf.pos.gz"
echo ""
echo "Next steps:"
echo "  - Run p06d (theta calculation)"
echo "  - Run p06c2 (beagle/MAF/IBS for PCA)"
echo "=================================================================="