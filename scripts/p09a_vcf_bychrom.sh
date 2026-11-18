#!/bin/bash
#SBATCH --job-name=vcf_chr
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=16:00:00
#SBATCH --output=logs/vcf_chr_%A_%a.log

# =============================================================================
# p09a_vcf_bychrom.sh — per-chromosome VCF/BCF generation for lcMLkin
# relatedness estimation. All 180 individuals, restricted to the existing
# LD-pruned site panel (12,922 SNPs from p06b) — the standard choice for
# pairwise relatedness (unlinked markers; genome-wide density adds nothing
# useful and multiplies compute for 16,110 pairwise comparisons).
#
# Same corrected filters as the rest of this session's pipeline (RagTag
# reference, -baq 1 -C 50 -uniqueOnly 1). Per-chromosome array — same
# no-checkpoint lesson as HWE: a single genome-wide ANGSD call has no resume,
# so any timeout costs the whole run; splitting means only one chromosome
# needs redoing.
#
# Modes:
#   --array-task  (ARRAY TASK) generate one chromosome's BCF
#   --concat      concatenate all 13 into one genome-wide VCF for lcMLkin
# =============================================================================

if [[ "${1:-}" == "--array-task" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    source ~/.bash_profile
    conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
    set -eu
    THREADS=${SLURM_CPUS_PER_TASK:-6}
    OUT_DIR=${OUT_DIR:-physalia/angsd_results}
    REFERENCE=${REFERENCE:-/mnt/parscratch/users/bi4og/ncbi_submission/plongirostris.fasta}
    CHROM_FILE=${CHROM_FILE:-$OUT_DIR/chromosomes.txt}
    SITES_FILE=${SITES_FILE:-$OUT_DIR/LDpruned_snps.list}
    BAM_LIST=${BAM_LIST:-physalia/sample_lists/bam_list_dedup_overlapclipped.list}
    MIND=${MIND:-83   # 50% of 166

    [[ -f "${SITES_FILE}.idx" ]] || { echo "ERROR: $SITES_FILE not indexed. Run: angsd sites index $SITES_FILE"; exit 1; }

    CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_FILE")
    [[ -n "$CHR" ]] || { echo "ERROR: no chromosome on line $SLURM_ARRAY_TASK_ID"; exit 1; }

    mkdir -p "$OUT_DIR/lcmlkin/vcf_parts"
    OUT_PREFIX="$OUT_DIR/lcmlkin/vcf_parts/ALL.${CHR}"

    if [[ -s "${OUT_PREFIX}.bcf" ]]; then
        if bcftools view "${OUT_PREFIX}.bcf" -h >/dev/null 2>&1; then
            echo "Already done: ${OUT_PREFIX}.bcf"; exit 0
        else
            echo "Existing BCF looks corrupt, regenerating"; rm -f "${OUT_PREFIX}.bcf"
        fi
    fi

    echo "=================================================================="
    echo "VCF/BCF  chr=$CHR  minInd=$MIND  threads=$THREADS"
    echo "=================================================================="

    angsd -bam "$BAM_LIST" -ref "$REFERENCE" \
        -r "$CHR" \
        -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 \
        -SNP_pval 1e-6 \
        -doBcf 1 -doPost 1 -doGeno 1 -doCounts 1 --ignore-RG 1 \
        -minMapQ 20 -minQ 20 -baq 1 -C 50 -uniqueOnly 1 \
        -minInd "$MIND" -setMinDepthInd 1 -remove_bads 1 -trim 0 \
        -sites "$SITES_FILE" \
        -P "$THREADS" \
        -out "$OUT_PREFIX" \
        2> "${OUT_PREFIX}.angsd.log"

    [[ -s "${OUT_PREFIX}.bcf" ]] || { echo "ERROR: BCF not created for $CHR"; tail -30 "${OUT_PREFIX}.angsd.log"; exit 1; }
    bcftools index -f "${OUT_PREFIX}.bcf"
    echo "Done: ${OUT_PREFIX}.bcf"
    exit 0
fi

if [[ "${1:-}" == "--concat" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    source ~/.bash_profile
    conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
    set -eu
    OUT_DIR=${OUT_DIR:-physalia/angsd_results}
    CHROM_FILE=${CHROM_FILE:-$OUT_DIR/chromosomes.txt}
    OUT_VCF="$OUT_DIR/lcmlkin/ALL.vcf.gz"

    PARTS=()
    while read -r c; do
        f="$OUT_DIR/lcmlkin/vcf_parts/ALL.${c}.bcf"
        [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
        PARTS+=("$f")
    done < "$CHROM_FILE"
    echo "Concatenating ${#PARTS[@]} chromosome BCFs..."

    bcftools concat "${PARTS[@]}" -O z -o "$OUT_VCF"
    bcftools index -ft "$OUT_VCF"

    N_SITES=$(bcftools view -H "$OUT_VCF" | wc -l)
    N_SAMPLES=$(bcftools query -l "$OUT_VCF" | wc -l)
    echo "Merged VCF: $N_SITES sites, $N_SAMPLES samples -> $OUT_VCF"
    [[ "$N_SITES" -gt 0 && "$N_SAMPLES" -eq 166 ]] || { echo "ERROR: sanity check failed (sites=$N_SITES samples=$N_SAMPLES, expected 166)"; exit 1; }
    exit 0
fi

echo "Usage: $0 --array-task|--concat [OUT_DIR=...] [...]"
exit 1
