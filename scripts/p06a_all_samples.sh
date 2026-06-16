#!/bin/bash
#SBATCH --job-name=p06a_all_samples
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p06a_all_samples_%A.log

# Load environment
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -euo pipefail

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

OUT_PREFIX="$OUT_DIR/all_samples"

mkdir -p "$OUT_DIR/logs"

# Validate inputs
if [[ ! -f "$BAM_LIST" ]]; then
    echo "ERROR: BAM list not found: $BAM_LIST"
    exit 1
fi

if [[ ! -f "$REFERENCE.fai" ]]; then
    echo "ERROR: Reference FAI not found at $REFERENCE.fai"
    exit 1
fi

N_SAMPLES=$(wc -l < "$BAM_LIST")
MIN_IND=$(echo "$N_SAMPLES * $MIN_IND_RATIO" | bc | awk '{print int($1+0.5)}')
[[ $MIN_IND -lt 1 ]] && MIN_IND=1
[[ $MIN_IND -gt $N_SAMPLES ]] && MIN_IND=$N_SAMPLES

echo "=================================================================="
echo "ANGSD: ALL Samples (Full Genome for LD Pruning)"
echo "=================================================================="
echo "Samples: $N_SAMPLES"
echo "Min individuals: $MIN_IND"
echo "Threads: $THREADS"
echo "=================================================================="

# Ancestral - use reference if not provided
if [[ -n "$ANCESTRAL" && -f "$ANCESTRAL" ]]; then
    ANC_OPT="-anc $ANCESTRAL"
    echo "Using ancestral genome: $ANCESTRAL"
else
    ANC_OPT="-anc $REFERENCE"
    echo "Using reference as ancestral (folded): $REFERENCE"
fi

# Run ANGSD genome-wide with STRICT filters for LD detection
echo ""
echo "Running ANGSD (genome-wide, strict filters for LD estimation)..."

if ! $ANGSD_BIN -bam "$BAM_LIST" -ref "$REFERENCE" $ANC_OPT \
    -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doCounts 1 \
    -minMapQ $MIN_MAPQ -minQ $MIN_Q -minInd $MIN_IND \
    -SNP_pval 1e-10 \
    -minMaf 0.05 \
    -setMinDepthInd 1 -nQueueSize 50 -P $THREADS \
    -out "$OUT_PREFIX" 2>&1 | tee "$OUT_PREFIX.angsd.log"; then
    echo "ERROR: ANGSD failed"
    tail -50 "$OUT_PREFIX.angsd.log"
    exit 1
fi

# Verify outputs
if [[ ! -f "$OUT_PREFIX.beagle.gz" ]]; then
    echo "ERROR: Beagle file not created"
    exit 1
fi

if [[ ! -f "$OUT_PREFIX.mafs.gz" ]]; then
    echo "ERROR: MAF file not created"
    exit 1
fi

echo "✓ ANGSD completed successfully"

# Show summary
echo ""
echo "Sites summary:"
grep -E "sites|Total number" "$OUT_PREFIX.angsd.log" | tail -10

echo ""
echo "=================================================================="
echo "ALL Samples ANGSD Complete"
echo "=================================================================="
echo "Output files:"
echo "  Beagle:  $OUT_PREFIX.beagle.gz"
echo "  MAF:     $OUT_PREFIX.mafs.gz"
echo "=================================================================="