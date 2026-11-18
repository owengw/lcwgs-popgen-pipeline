#!/bin/bash
#SBATCH --job-name=pair_het
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=1-44%22
#SBATCH --output=logs/p07t_pair_het_%A_%a.log

# =============================================================================
# p07t_pair_het.sh - batch test for individual heterozygosity.
# Same method as p07a_individual_het.sh (1-sample folded SFS), run on the 11
# animals sequenced in BOTH batches, on chr1 only (to be quick), twice:
#   tasks  1-22: -setMinDepthInd 1  (the setting used in p07a)
#   tasks 23-44: -setMinDepthInd 5  (only sites with enough reads to see a het)
# If the two copies of each animal agree at one setting, that setting is
# robust to the batch difference. Compare with p07t_compare.R.
# =============================================================================

_SAVED_ARGS=("$@")
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -- "${_SAVED_ARGS[@]}"
set -eu

OUT_DIR="physalia/angsd_results"
REFERENCE="/mnt/parscratch/users/bi4og/ncbi_submission/plongirostris.fasta"
BAM_LIST="physalia/sample_lists/bam_list_batch_pairs.list"
REGION="chr1_RagTag_pilon"
THREADS=${SLURM_CPUS_PER_TASK:-2}

N=$(wc -l < "$BAM_LIST")
[[ "$N" -gt 0 ]] || { echo "ERROR: $BAM_LIST is empty - run p07t_make_pair_list.sh"; exit 1; }
T=${SLURM_ARRAY_TASK_ID}
if [[ $T -le $N ]]; then MIND=1; LINE=$T; else MIND=5; LINE=$((T - N)); fi
[[ $LINE -le $N ]] || { echo "Task $T beyond 2 x $N BAMs - nothing to do"; exit 0; }

BAM=$(sed -n "${LINE}p" "$BAM_LIST")
SAMPLE_ID=$(basename "$BAM" | sed 's/_merged.*//')
HET_DIR="$OUT_DIR/het_batch_test/mindepth${MIND}"
mkdir -p "$HET_DIR" logs
SAF_PREFIX="$HET_DIR/${SAMPLE_ID}"
HET_FILE="${SAF_PREFIX}.het"
[[ -s "$HET_FILE" ]] && { echo "Already done: $HET_FILE"; cat "$HET_FILE"; exit 0; }
echo "Sample $SAMPLE_ID  setMinDepthInd=$MIND  region=$REGION"

angsd -i "$BAM" -ref "$REFERENCE" -anc "$REFERENCE" -r "$REGION" \
    -GL 2 -doSaf 1 -doCounts 1 \
    -minMapQ 20 -minQ 20 -baq 1 -C 50 -uniqueOnly 1 \
    -setMinDepthInd $MIND -setMaxDepthInd 50 \
    -remove_bads 1 -trim 0 -P $THREADS -nQueueSize 50 \
    -out "$SAF_PREFIX" 2> "${SAF_PREFIX}.angsd.log"
[[ -s "${SAF_PREFIX}.saf.idx" ]] || { echo "ERROR: SAF not created"; tail -30 "${SAF_PREFIX}.angsd.log"; exit 1; }

realSFS "${SAF_PREFIX}.saf.idx" -fold 1 -P $THREADS > "${SAF_PREFIX}.sfs" 2> "${SAF_PREFIX}.sfs.log"
awk -v s="$SAMPLE_ID" -v m="$MIND" 'BEGIN{OFS="\t"; print "sample_id","min_depth","heterozygosity","n_sites_total"}
     {t=$1+$2+$3; if (t>0) print s, m, sprintf("%.8f",$2/t), sprintf("%.1f",t)}' "${SAF_PREFIX}.sfs" > "$HET_FILE"
cat "$HET_FILE"
rm -f "${SAF_PREFIX}.saf.gz" "${SAF_PREFIX}.saf.pos.gz"   # SAFs not needed afterwards
