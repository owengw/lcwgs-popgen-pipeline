#!/bin/bash
#SBATCH --job-name=p06c2_pca
#SBATCH --cpus-per-task=40
#SBATCH --mem=160G
#SBATCH --time=96:00:00
#SBATCH --output=logs/p06c2_pca_%A.log

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
THREADS=${THREADS:-40}
MIN_MAPQ=${MIN_MAPQ:-10}
MIN_Q=${MIN_Q:-13}
MIN_IND_RATIO=${MIN_IND_RATIO:-0.20}

# Validate inputs
if [[ ! -f "$BAM_LIST" ]]; then
    echo "ERROR: BAM list not found: $BAM_LIST"
    exit 1
fi

if [[ ! -f "$REFERENCE.fai" ]]; then
    echo "ERROR: Reference index not found: $REFERENCE.fai"
    exit 1
fi

# Validate LD-pruned files exist
if [[ ! -f "$OUT_DIR/LDpruned_snps.list" ]]; then
    echo "ERROR: LD-pruned sites file not found: $OUT_DIR/LDpruned_snps.list"
    echo "Run step 1 (LD pruning) first"
    exit 1
fi

if [[ ! -f "$OUT_DIR/LDpruned_contigs.txt" ]]; then
    echo "ERROR: LD-pruned contigs file not found: $OUT_DIR/LDpruned_contigs.txt"
    echo "Run step 1 (LD pruning) first"
    exit 1
fi

# Extract population name
POP_NAME=$(basename "$BAM_LIST" .list)
OUT_PREFIX="$OUT_DIR/${POP_NAME}"

mkdir -p "$OUT_DIR/logs"

# Calculate minInd
N_SAMPLES=$(wc -l < "$BAM_LIST")
MIN_IND=$(echo "$N_SAMPLES * $MIN_IND_RATIO" | bc | awk '{print int($1+0.5)}')
[[ $MIN_IND -lt 1 ]] && MIN_IND=1
[[ $MIN_IND -gt $N_SAMPLES ]] && MIN_IND=$N_SAMPLES

echo "=================================================================="
echo "ANGSD: $POP_NAME (PCA: Beagle/MAF/IBS)"
echo "=================================================================="
echo "Population: $POP_NAME"
echo "Samples: $N_SAMPLES"
echo "Min individuals: $MIN_IND"
echo "Threads: $THREADS"
echo "Using LD-pruned sites for PCA (avoids LD inflation)"
echo "=================================================================="

# Ancestral genome handling
if [[ -n "$ANCESTRAL" && -f "$ANCESTRAL" ]]; then
    ANC_OPT="-anc $ANCESTRAL"
    echo "Using ancestral genome: $ANCESTRAL"
else
    ANC_OPT="-anc $REFERENCE"
    echo "Using reference as ancestral: $REFERENCE"
fi

# -------------------------------
# ANGSD Run: Beagle, MAF, IBS with LD-pruned sites
# Using -doMajorMinor 3 (from sites file) for memory efficiency
# Using -GL 1 (SAMtools) instead of -GL 2 (GATK) for lighter processing
# -------------------------------
echo ""
echo "Running ANGSD (Beagle/MAF/IBS for PCA)..."
echo "Using -doMajorMinor 3: major/minor from LD-pruned sites file"
echo "Using -GL 1: SAMtools likelihood model"

if ! $ANGSD_BIN -bam "$BAM_LIST" -ref "$REFERENCE" $ANC_OPT \
    -GL 1 -doGlf 2 -doMajorMinor 3 -doMaf 1 -doPost 1 -doCounts 1 \
    -doIBS 1 -makeMatrix 1 -doCov 1 \
    -minMapQ $MIN_MAPQ -minQ $MIN_Q -minInd $MIN_IND \
    -nQueueSize 50 -P $THREADS \
    -sites "$OUT_DIR/LDpruned_snps.list" \
    -rf "$OUT_DIR/LDpruned_contigs.txt" \
    -out "$OUT_PREFIX" 2>&1 | tee "$OUT_PREFIX.angsd_pca.log"; then
    echo "ERROR: ANGSD PCA analysis failed"
    tail -50 "$OUT_PREFIX.angsd_pca.log"
    exit 1
fi

# Verify outputs
MISSING_FILES=0
if [[ ! -f "$OUT_PREFIX.beagle.gz" ]]; then
    echo "WARNING: Beagle file not created"
    MISSING_FILES=$((MISSING_FILES + 1))
fi

if [[ ! -f "$OUT_PREFIX.mafs.gz" ]]; then
    echo "WARNING: MAF file not created"
    MISSING_FILES=$((MISSING_FILES + 1))
fi

if [[ ! -f "$OUT_PREFIX.covMat" ]]; then
    echo "WARNING: Covariance matrix not created"
    MISSING_FILES=$((MISSING_FILES + 1))
fi

if [[ ! -f "$OUT_PREFIX.ibsMat" ]]; then
    echo "WARNING: IBS matrix not created"
    MISSING_FILES=$((MISSING_FILES + 1))
fi

if [[ $MISSING_FILES -gt 0 ]]; then
    echo "ERROR: $MISSING_FILES output files missing"
    exit 1
fi

echo "✓ ANGSD PCA analysis completed successfully"

# Show summary
echo ""
echo "PCA analysis summary:"
grep -E "sites|Total number" "$OUT_PREFIX.angsd_pca.log" | tail -10

echo ""
echo "=================================================================="
echo "ANGSD PCA Complete: $POP_NAME"
echo "=================================================================="
echo "Output files:"
echo "  Beagle:     $OUT_PREFIX.beagle.gz"
echo "  MAF:        $OUT_PREFIX.mafs.gz"
echo "  Cov matrix: $OUT_PREFIX.covMat"
echo "  IBS matrix: $OUT_PREFIX.ibsMat"
echo ""
echo "Use covMat or ibsMat for PCA visualization"
echo "Example PCA in R:"
echo "  cov <- as.matrix(read.table('$OUT_PREFIX.covMat'))"
echo "  pca <- eigen(cov)"
echo "  plot(pca\$vectors[,1], pca\$vectors[,2])"
echo "=================================================================="