#!/usr/bin/env bash
#SBATCH --job-name=angsd_vcf
#SBATCH --output=/mnt/parscratch/users/bi4og/physalia/dsuite/logs/angsd_vcf_%j.log
#SBATCH --mem=250G
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8

set -euo pipefail

# =============================================================================
# Generate VCF from ANGSD for Dsuite ABBA-BABA analysis
# All 7 populations: CAI, NS, SB, SI, COI, DH, SP
# Sites: 10M random site list
# =============================================================================

BASE="/mnt/parscratch/users/bi4og/physalia"
ANGSD_DIR="${BASE}/angsd_results"
GENOME="${BASE}/../genome/Plong_genome_flye.fasta"
SITES="${ANGSD_DIR}/random_10M_sites.list"
CONTIGS="${ANGSD_DIR}/random_10M_contigs.txt"
OUT_DIR="${BASE}/dsuite"
OUT_PREFIX="${OUT_DIR}/pllong_all_pops"

mkdir -p "${OUT_DIR}/logs"

if [[ -f "${OUT_PREFIX}.bcf" ]]; then
    echo "SKIPPING ANGSD — BCF already exists: ${OUT_PREFIX}.bcf"
else

echo "ANGSD VCF generation for Dsuite"
echo "Started: $(date)"

# Activate conda
source /mnt/community/Genomics/apps/miniforge/miniforge3/etc/profile.d/conda.sh
set +u
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -u

if ! command -v angsd &>/dev/null; then
    echo "ERROR: angsd not found"
    exit 1
fi
echo "ANGSD: $(angsd 2>&1 | grep 'angsd version' | head -1)"

# =============================================================================
# Build combined BAM list (all 7 populations)
# Duplicates already removed by the DUP_PATTERN in previous analyses
# Use pop_map_ALL.list which contains all 198 samples (pre-dedup)
# For VCF we use the per-population lists to exclude duplicates
# =============================================================================

POPS=(CAI NS SB SI COI DH SP)
POP_LIST_DIR="${ANGSD_DIR}/pop_map_lists"
COMBINED_BAMS="${OUT_DIR}/all_pops_bams.list"

> "${COMBINED_BAMS}"
for POP in "${POPS[@]}"; do
    BAM_LIST="${POP_LIST_DIR}/pop_map_${POP}.list"
    if [[ ! -f "${BAM_LIST}" ]]; then
        echo "ERROR: BAM list not found: ${BAM_LIST}"
        exit 1
    fi
    N=$(wc -l < "${BAM_LIST}")
    echo "  ${POP}: ${N} individuals"
    cat "${BAM_LIST}" >> "${COMBINED_BAMS}"
done

N_TOTAL=$(wc -l < "${COMBINED_BAMS}")
echo "Total individuals: ${N_TOTAL}"

# =============================================================================
# Run ANGSD to generate VCF
# -doGeno -1: output genotype likelihoods
# -doVcf 1: output VCF format
# -doPost 1: compute posterior genotype probabilities using HWE prior
# -doMajorMinor 1: infer major/minor from GL
# -doMaf 1: compute allele frequencies
# -SNP_pval 1e-6: only output SNPs
# -minInd: require data in at least half of all individuals
# =============================================================================

MIN_IND=99  # half of 198 actual samples
echo "minInd: ${MIN_IND}"

angsd \
    -bam          "${COMBINED_BAMS}" \
    -ref          "${GENOME}" \
    -GL           1 \
    -doGlf        2 \
    -doGeno       1 \
    -doBcf        1 \
    --ignore-RG   1 \
    -doPost       1 \
    -doMajorMinor 1 \
    -doMaf        1 \
    -doCounts     1 \
    -SNP_pval     1e-6 \
    -minMapQ      10 \
    -minQ         20 \
    -minInd       "${MIN_IND}" \
    -setMinDepthInd 1 \
    -P            4 \
    -nQueueSize   1 \
    -remove_bads  1 \
    -trim         0 \
    -sites        "${SITES}" \
    -rf           "${CONTIGS}" \
    -out          "${OUT_PREFIX}"

if [[ ! -f "${OUT_PREFIX}.vcf.gz" ]] && [[ ! -f "${OUT_PREFIX}.vcf" ]]; then
    echo "ERROR: VCF not produced"
    exit 1
fi

fi
# Compress and index if not already done
# Convert BCF to bgzipped VCF and index
echo "Converting BCF to VCF..."
LD_PRELOAD="" bcftools view "${OUT_PREFIX}.bcf" | bgzip > "${OUT_PREFIX}.vcf.gz"
tabix -p vcf "${OUT_PREFIX}.vcf.gz"

N_SNPS=$(zcat "${OUT_PREFIX}.vcf.gz" 2>/dev/null | grep -v "^#" | wc -l || \
         grep -v "^#" "${OUT_PREFIX}.vcf" | wc -l)

echo ""
echo "VCF generation complete: $(date)"
echo "SNPs in VCF: ${N_SNPS}"
ls -lh "${OUT_PREFIX}".vcf* 2>/dev/null

# =============================================================================
# Generate Dsuite population assignment file
# Format: SampleID <tab> Population
# Sample IDs must match VCF sample names exactly
# =============================================================================

POP_FILE="${OUT_DIR}/dsuite_populations.txt"
> "${POP_FILE}"

for POP in "${POPS[@]}"; do
    BAM_LIST="${POP_LIST_DIR}/pop_map_${POP}.list"
    # Use DH as outgroup since it is the most differentiated single individual
    POP_LABEL="${POP}"
    [[ "${POP}" == "DH" ]] && POP_LABEL="Outgroup"
    while IFS= read -r bam; do
        echo -e "${bam}\t${POP_LABEL}" >> "${POP_FILE}"
    done < "${BAM_LIST}"
done

echo "Population file written: ${POP_FILE}"
echo "  Lines: $(wc -l < ${POP_FILE})"
head -5 "${POP_FILE}"

echo "All done: $(date)"
