#!/bin/bash
#SBATCH --job-name=angsd_vcf
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-4
#SBATCH --output=logs/angsd_vcf_%A_%a.log

# =============================================================================
# Generate per-population VCFs for lcMLkin relatedness analysis
#
# Uses ANGSD with hard genotype calling (-dogeno 1, -doPost 1)
# Filters: minMapQ 20, minQ 20, minInd 50% of population, SNP p<1e-6
# Output: per-population VCF with hard genotype calls
#
# Run BEFORE p09_lcmlkin.sh
# =============================================================================

source ~/.bash_profile
module load ANGSD/0.940-GCC-12.2.0 2>/dev/null || \
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -euo pipefail
set -x

# =============================================================================
# Variables
# =============================================================================
OUT_DIR=${OUT_DIR:-physalia/angsd_results}
REF="genome/Plong_genome_flye.fasta"
POPULATIONS=(CAI NS SB SI)
POP=${POPULATIONS[$((SLURM_ARRAY_TASK_ID - 1))]}
POP_LIST="$OUT_DIR/pop_map_lists/pop_map_${POP}.list"
VCF_OUT="$OUT_DIR/lcmlkin/${POP}/${POP}_angsd"
THREADS=$SLURM_CPUS_PER_TASK
EXCLUDE="o_merged\|T_merged\|SI43\|SI45\|SI83"

mkdir -p "$OUT_DIR/lcmlkin/$POP"
mkdir -p logs

echo "=========================================="
echo "Generating VCF: $POP"
echo "Array task: $SLURM_ARRAY_TASK_ID / 4"
echo "=========================================="

# Validate inputs
[[ ! -f "$REF" ]]      && echo "ERROR: Reference not found: $REF"       && exit 1
[[ ! -f "$POP_LIST" ]] && echo "ERROR: BAM list not found: $POP_LIST"   && exit 1

# Build filtered BAM list (exclude duplicates/failures)
FILTERED_LIST="$OUT_DIR/lcmlkin/$POP/${POP}_bams_filtered.list"
grep -v "$EXCLUDE" "$POP_LIST" > "$FILTERED_LIST"
N_BAMS=$(wc -l < "$FILTERED_LIST")
MIN_IND=$(echo "$N_BAMS * 0.5" | bc | awk '{print int($1)}')
echo "  BAMs: $N_BAMS, minInd: $MIN_IND"

# =============================================================================
# Run ANGSD
# =============================================================================
angsd \
    -bam "$FILTERED_LIST" \
    -ref "$REF" \
    -out "$VCF_OUT" \
    -nThreads "$THREADS" \
    -remove_bads 1 \
    -C 50 \
    -minMapQ 20 \
    -minQ 20 \
    -minInd "$MIN_IND" \
    -GL 2 \
    -dogeno 1 \
    -doPost 1 \
    -doMaf 1 \
    -doMajorMinor 1 \
    -SNP_pval 1e-6 \
    --ignore-RG 1 \
    -docounts 1 \
    -doBcf 1 \
    2> "$OUT_DIR/lcmlkin/$POP/${POP}_angsd.log"

# Check output
if [[ ! -f "${VCF_OUT}.bcf" ]]; then
    echo "ERROR: ANGSD VCF generation failed"
    tail -10 "$OUT_DIR/lcmlkin/$POP/${POP}_angsd.log"
    exit 1
fi

# Convert BCF to VCF.gz
module load BCFtools/1.17-GCC-12.2.0 2>/dev/null || \
    BCFTOOLS=/mnt/parscratch/users/bi4og/conda_envs/owenspopgen/bin/bcftools

${BCFTOOLS:-bcftools} view \
    --output-type z \
    --output "${VCF_OUT}.vcf.gz" \
    "${VCF_OUT}.bcf"
${BCFTOOLS:-bcftools} index -t "${VCF_OUT}.vcf.gz"

N_SITES=$(${BCFTOOLS:-bcftools} stats "${VCF_OUT}.vcf.gz" | \
          grep "^SN.*number of SNPs" | awk '{print $NF}')
N_SAMP=$(${BCFTOOLS:-bcftools} query -l "${VCF_OUT}.vcf.gz" | wc -l)

echo "  VCF: $N_SAMP samples, $N_SITES SNPs"
echo "  Output: ${VCF_OUT}.vcf.gz"

echo "=========================================="
echo "ANGSD VCF complete: $POP"
echo "=========================================="
