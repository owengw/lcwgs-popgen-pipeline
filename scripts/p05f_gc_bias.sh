#!/bin/bash
#SBATCH --job-name=gc_bias
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=1-20%14
#SBATCH --output=logs/gc_bias_%A_%a.log

# =============================================================================
# p05f_gc_bias.sh — GC-content bias check via samtools stats (GCF/GCL),
# replacing the earlier Picard CollectGcBiasMetrics approach. samtools stats
# reads only the BAM itself, with no separate reference-dictionary matching
# requirement — sidesteps the unplaced-contig dictionary drift entirely.
#
# GCF/GCL report the % of first/second-in-pair reads at each GC-content bin
# (0-100), directly from the reads — a standard measure of compositional
# bias, appropriate for detecting WGA-driven GC skew.
# =============================================================================

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -eu

OUT_DIR="physalia/angsd_results/amp_bias/gc_bias"
mkdir -p "$OUT_DIR"

PAIRS=(
  "289o-NS40o_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:283-NS40_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:NS40"
  "294o-SI7o_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:221-SI7_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:SI7"
  "296o-NS14o_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:227-NS14_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:NS14"
  "297o-NS41o_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:291-NS41_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:NS41"
  "315o-SI18o_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:198-SI18_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:SI18"
  "332o-CAI6o_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:327-CAI6_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:CAI6"
  "337o-SI59o_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:233-SI59_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:SI59"
  "351o-NS1o_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:258-NS1_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:NS1"
  "352o-NS30o_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:236-NS30_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:NS30"
  "359o-NS32o_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:244-NS32_merged_pe_bt2_plongirostris_minq20_sorted_dedup_overlapclipped.bam:NS32"
)

IDX=$SLURM_ARRAY_TASK_ID
# field 1 = original unamplified DNA ("o"), field 2 = amplified library of the same extract
if [[ $IDX -le 10 ]]; then
    PAIR="${PAIRS[$((IDX-1))]}"; BAM=$(echo "$PAIR" | cut -d: -f1); IND=$(echo "$PAIR" | cut -d: -f3)
    LABEL="${IND}_non_amplified"
else
    PAIR="${PAIRS[$((IDX-11))]}"; BAM=$(echo "$PAIR" | cut -d: -f2); IND=$(echo "$PAIR" | cut -d: -f3)
    LABEL="${IND}_amplified"
fi

if [[ -f "physalia/merged/$BAM" ]]; then
    BAM_PATH="physalia/merged/$BAM"
elif [[ -f "physalia/merged_excluded/$BAM" ]]; then
    BAM_PATH="physalia/merged_excluded/$BAM"
else
    echo "ERROR: BAM not found in merged/ or merged_excluded/: $BAM"; exit 1
fi

echo "=================================================================="
echo "GC bias (samtools stats): $LABEL  ($BAM_PATH)"
echo "=================================================================="

OUT="$OUT_DIR/${LABEL}.samtools_stats.txt"
if [[ -s "$OUT" ]]; then
    echo "Already done: $OUT"; exit 0
fi

samtools stats "$BAM_PATH" > "$OUT"

# quick sanity check
N_GCF=$(grep -c "^GCF" "$OUT")
[[ "$N_GCF" -gt 0 ]] || { echo "ERROR: no GCF lines in output"; exit 1; }
echo "Done: $LABEL ($N_GCF GCF bins)"
