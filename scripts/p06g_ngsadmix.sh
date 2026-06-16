#!/bin/bash
#SBATCH --job-name=ngsadmix
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-200
#SBATCH --output=logs/ngsadmix_%A_%a.log

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -euo pipefail
set -x

OUT_DIR=${OUT_DIR:-physalia/angsd_results}
BEAGLE="$OUT_DIR/pcangsd_ALL.beagle.gz"
NGSADMIX_DIR="$OUT_DIR/ngsadmix"
THREADS=4

mkdir -p "$NGSADMIX_DIR"

# Calculate K and rep from array task ID
K=$(( (SLURM_ARRAY_TASK_ID - 1) / 20 + 1 ))
REP=$(( (SLURM_ARRAY_TASK_ID - 1) % 20 + 1 ))

SEED=$((K * 1000 + REP))
OUT_PREFIX="$NGSADMIX_DIR/K${K}_rep${REP}/ngsadmix_K${K}_rep${REP}"

mkdir -p "$NGSADMIX_DIR/K${K}_rep${REP}"

echo "=========================================="
echo "NGSadmix: K=$K Rep=$REP Seed=$SEED"
echo "=========================================="

NGSadmix \
    -likes "$BEAGLE" \
    -K $K \
    -seed $SEED \
    -o "$OUT_PREFIX" \
    -P $THREADS \
    -minMaf 0.05

# Check output and extract log likelihood
if [[ -f "${OUT_PREFIX}.qopt" ]]; then
    LOGLIKE=$(grep "best like=" "${OUT_PREFIX}.log" | grep -oP 'best like=\K[-0-9.]+')
    echo "K=$K Rep=$REP loglike=$LOGLIKE" >> "$NGSADMIX_DIR/loglikelihoods.txt"
    echo "NGSadmix completed: K=$K Rep=$REP loglike=$LOGLIKE"
else
    echo "ERROR: NGSadmix failed for K=$K Rep=$REP"
    exit 1
fi