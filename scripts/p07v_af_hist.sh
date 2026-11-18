#!/bin/bash
#SBATCH --job-name=af_hist
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/p07v_af_hist_%A_%a.log

# =============================================================================
# p07v_af_hist.sh - what does the second allele look like at well-covered sites?
# For each BAM in the p07u list: at every site with 10-30 reads on the first 20 Mb
# of chr1, the fraction of reads carrying the commonest non-reference allele.
# The shape tells the cause apart:
#   clean single animal : a peak at 0.5 (hets) and almost nothing at 0.05-0.3
#   DNA from 2 animals  : extra peaks at about c/2 and c, where c is the share of
#                         the second animal's DNA (e.g. 0.1 and 0.2 for 20%)
#   allele dropout /    : hets smeared from 0.1 to 0.4 with no clear peak
#   amplification bias
# Submit from /mnt/parscratch/users/bi4og (uses the list made by p07u --list):
#   sbatch --array=1-N scripts/p07v_af_hist.sh       (N = lines in the p07u list)
#   Rscript scripts/p07v_af_hist.R                    (plot + table, after it finishes)
# =============================================================================
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -eu

OUT_DIR="physalia/angsd_results/run_qc/af_hist"
LIST="physalia/sample_lists/bam_list_run_qc.list"
REFERENCE="/mnt/parscratch/users/bi4og/ncbi_submission/plongirostris.fasta"
REGION="chr1_RagTag_pilon:1-20000000"
mkdir -p "$OUT_DIR" logs

BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[[ -n "$BAM" ]] || { echo "no BAM on line ${SLURM_ARRAY_TASK_ID} of $LIST"; exit 1; }
ID=$(basename "$BAM" | sed 's/_merged.*//')
echo "$ID  $REGION"

# AD = reads per allele (ref first). Keep the commonest alternative allele if it has >= 2 reads; bins of 0.025.
bcftools mpileup -f "$REFERENCE" -r "$REGION" -q 20 -Q 20 -a AD,DP -Ou "$BAM" 2>/dev/null \
  | bcftools query -i 'FMT/DP>=10 && FMT/DP<=30' -f '[%AD]\n' \
  | awk -F, -v id="$ID" '
      { t=0; mx=0; for(i=1;i<=NF;i++){ t+=$i; if(i>1 && $i>mx) mx=$i }
        n++; if (t==0 || mx==0) { zero++; next }
        if (mx==1) { single++; next }            # one read only: usually a sequencing error
        b=int((mx/t)/0.025); if (b>39) b=39; h[b]++ }
      END { print "sample\tbin_low\tsites"; print id "\t-1\t" zero+0; print id "\t-2\t" single+0
            for(b=0;b<40;b++) printf "%s\t%.3f\t%d\n", id, b*0.025, h[b]+0 }' > "$OUT_DIR/$ID.afhist.tsv"
echo "done: $OUT_DIR/$ID.afhist.tsv"
