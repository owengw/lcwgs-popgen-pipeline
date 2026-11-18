#!/bin/bash
#SBATCH --job-name=pca_beagle_chr
#SBATCH --cpus-per-task=6
#SBATCH --mem=128G
#SBATCH --time=30:00:00
#SBATCH --output=logs/pca_beagle_chr_%A_%a.log

# =============================================================================
# p06h_pca_beagle_bychrom.sh — REGENERATES pcangsd_ALL.beagle.gz with the
# BAQ/C/uniqueOnly correction. The original (Jul 29) predates the fix
# entirely — confirmed via .arg: no -baq, no -C, no -uniqueOnly. This beagle
# feeds PCA, PCAngsd admixture (K=2-4), and PCAngsd individual F (used both
# for inbreeding results AND as cross-validation for the sex-determination
# chr5 finding) — all need regenerating from this corrected version.
#
# Same LD-pruned site panel (LDpruned_snps.list, -doMajorMinor 3 since the
# file carries major/minor alleles), 180 samples jointly (heaviest joint
# operation after the full lcMLkin VCF — pre-set to 128G/30h based on that
# scaling, not incremental guessing).
#
# Modes: --array-task (one chromosome), --concat (merge all 13)
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
    BAM_LIST=${BAM_LIST:-physalia/angsd_results/pop_map_lists/pop_map_ALL.list}

    CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHROM_FILE")
    [[ -n "$CHR" ]] || { echo "ERROR: no chromosome on line $SLURM_ARRAY_TASK_ID"; exit 1; }

    mkdir -p "$OUT_DIR/beagle_parts_corrected"
    OUT_PREFIX="$OUT_DIR/beagle_parts_corrected/pcangsd_ALL.${CHR}"

    if [[ -s "${OUT_PREFIX}.beagle.gz" ]]; then
        gzip -t "${OUT_PREFIX}.beagle.gz" 2>/dev/null && { echo "Already done: ${OUT_PREFIX}.beagle.gz"; exit 0; }
        echo "Existing beagle truncated, regenerating"; rm -f "${OUT_PREFIX}.beagle.gz"
    fi

    echo "=================================================================="
    echo "PCA beagle (corrected)  chr=$CHR  threads=$THREADS"
    echo "=================================================================="

    angsd -bam "$BAM_LIST" -ref "$REFERENCE" \
        -r "$CHR" \
        -GL 2 -doGlf 2 -doMajorMinor 3 -doCounts 1 -doMaf 1 \
        -minMapQ 30 -minQ 20 -baq 1 -C 50 -uniqueOnly 1 \
        -sites "$SITES_FILE" \
        -nQueueSize 50 -P "$THREADS" \
        -out "$OUT_PREFIX" \
        2> "${OUT_PREFIX}.angsd.log"

    gzip -t "${OUT_PREFIX}.beagle.gz" || { echo "ERROR: beagle corrupt for $CHR"; exit 1; }
    echo "Done: ${OUT_PREFIX}.beagle.gz"
    exit 0
fi

if [[ "${1:-}" == "--concat" ]]; then
    shift
    for kv in "$@"; do export "$kv"; done
    source ~/.bash_profile
    set -eu
    OUT_DIR=${OUT_DIR:-physalia/angsd_results}
    CHROM_FILE=${CHROM_FILE:-$OUT_DIR/chromosomes.txt}
    OUT_BEAGLE="$OUT_DIR/pcangsd_ALL_corrected.beagle.gz"

    PARTS=()
    while read -r c; do
        f="$OUT_DIR/beagle_parts_corrected/pcangsd_ALL.${c}.beagle.gz"
        [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
        gzip -t "$f" || { echo "ERROR: corrupt $f"; exit 1; }
        PARTS+=("$f")
    done < "$CHROM_FILE"
    echo "All ${#PARTS[@]} chromosome parts present and valid"

    { zcat "${PARTS[0]}" | head -1
      for f in "${PARTS[@]}"; do zcat "$f" | tail -n +2; done
    } | gzip > "$OUT_BEAGLE"

    N_SITES=$(zcat "$OUT_BEAGLE" | tail -n +2 | wc -l)
    N_SAMP=$(zcat "$OUT_BEAGLE" | head -1 | awk '{print (NF-3)/3}')
    echo "Merged: $N_SITES sites, $N_SAMP individuals -> $OUT_BEAGLE"
    [[ "$N_SAMP" -eq 166 ]] || { echo "ERROR: expected 166 individuals, got $N_SAMP"; exit 1; }
    exit 0
fi

echo "Usage: $0 --array-task|--concat [OUT_DIR=...]"
exit 1
