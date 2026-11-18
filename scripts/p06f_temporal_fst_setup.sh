#!/bin/bash
#SBATCH --job-name=temporal_fst_setup
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/temporal_fst_setup_%j.log

# =============================================================================
# p06f_temporal_fst_setup.sh — builds year-split sample lists (2023 vs 2024)
# for ALL individuals and for CAI/NS/SB populations (the three with enough
# samples per year, per the methods). Then launch the SAME corrected
# per-chromosome FST pipeline (p06e_fst_bychrom.sh) on each year-pair,
# exactly as used for the spatial comparison — this reuses the fixed
# -nSites-chunked, BAQ/C-corrected SAF/2D-SFS machinery, not a rerun of the
# old pre-fix temporal FST.
# =============================================================================

source ~/.bash_profile
set -eu

OUT_DIR="${OUT_DIR:-physalia/angsd_results}"
POP_LIST_DIR="$OUT_DIR/pop_map_lists"
YEAR_DIR="$OUT_DIR/pop_map_lists_byyear"
METADATA="$OUT_DIR/metadata.tsv"
mkdir -p "$YEAR_DIR"

# metadata.tsv: sample_id  <tab> ... <tab> population <tab> year (columns per
# earlier p07c usage: sample_prefix, site, population, year — confirm order)
for GROUP in ALL CAI NS SB; do
    FULL_LIST="$POP_LIST_DIR/pop_map_${GROUP}.list"
    [[ -f "$FULL_LIST" ]] || { echo "WARNING: $FULL_LIST missing, skipping $GROUP"; continue; }

    for YEAR in 2023 2024; do
        OUT="$YEAR_DIR/pop_map_${GROUP}_${YEAR}.list"
        # match each BAM's sample prefix against metadata's year column
        > "$OUT"
        while read -r bam; do
            sample=$(basename "$bam" | sed 's/_merged.*//')
            y=$(awk -F'\t' -v s="$sample" '$1==s{print $NF}' "$METADATA")
            [[ "$y" == "$YEAR" ]] && echo "$bam" >> "$OUT"
        done < "$FULL_LIST"
        N=$(wc -l < "$OUT")
        echo "$GROUP $YEAR: $N samples -> $OUT"
    done
done
