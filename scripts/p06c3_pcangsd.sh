#!/bin/bash
#SBATCH --job-name=p06c_pcangsd
#SBATCH --cpus-per-task=40
#SBATCH --mem=160G
#SBATCH --time=96:00:00
#SBATCH --output=logs/%x_%A.log

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -euo pipefail
set -x

BAM_LIST=${BAM_LIST}
REFERENCE=${REFERENCE}
OUT_DIR=${OUT_DIR}
THREADS=${THREADS}

OUT_PREFIX="$OUT_DIR/pcangsd_ALL"

echo "Running ANGSD to generate Beagle file..."

angsd \
-bam "$BAM_LIST" \
-ref "$REFERENCE" \
-GL 2 \
-doGlf 2 \
-doMajorMinor 1 \
-sites "$OUT_DIR/LDpruned_snps.list" \
-rf "$OUT_DIR/LDpruned_contigs.txt" \
-doCounts 1 \
-doMaf 1 \
-minMapQ 30 \
-minQ 20 \
-nQueueSize 50 \
-P $THREADS \
-out "$OUT_PREFIX"

echo "Running PCAngsd..."

pcangsd.py \
-beagle "$OUT_PREFIX.beagle.gz" \
-o "$OUT_PREFIX" \
-threads $THREADS