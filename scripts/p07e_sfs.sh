#!/bin/bash
#SBATCH --job-name=sfs_analysis
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/p07e_sfs_%j.log

# =============================================================================
# p07e_sfs.sh — RagTag rewrite. Summarises singleton/doubleton proportions
# from the genome-wide summed SFS files already produced by
# p06d_theta_bychrom.sh --sfs-sum (pop_map_<POP>.sfs). No recomputation.
# Populations: CAI, NS, SB, SI (the four with full all-sites theta chains).
# Sample sizes read from pop_map_lists/, matching the current 180-sample set.
# =============================================================================

source ~/.bash_profile
set -eu

OUT_DIR=${OUT_DIR:-physalia/angsd_results}
POP_LIST_DIR="$OUT_DIR/pop_map_lists"
SFS_OUT="$OUT_DIR/sfs"
mkdir -p "$SFS_OUT" logs

POPULATIONS=(CAI NS SB SI)

SUMMARY="$SFS_OUT/sfs_summary.txt"
echo "SFS Summary Statistics (all-sites, genome-wide)" > "$SUMMARY"
echo "=================================================" >> "$SUMMARY"
echo "" >> "$SUMMARY"

for POP in "${POPULATIONS[@]}"; do
    SFS_FILE="$OUT_DIR/pop_map_${POP}.sfs"
    LIST_FILE="$POP_LIST_DIR/pop_map_${POP}.list"

    if [[ ! -s "$SFS_FILE" ]]; then
        echo "WARNING: SFS file not found for $POP ($SFS_FILE)"
        continue
    fi
    [[ -f "$LIST_FILE" ]] || { echo "WARNING: no BAM list for $POP"; continue; }

    N_IND=$(wc -l < "$LIST_FILE")
    echo "Processing $POP (n=$N_IND individuals, $((N_IND*2)) chromosomes)..."

    # folded SFS: entries 1..N_IND (index 0 = invariant, excluded)
    cat "$SFS_FILE" | tr ' ' '\n' | awk -v n="$N_IND" 'NF>0 && NR>1 && NR<=(n+1) {print}' \
        > "$SFS_OUT/${POP}_sfs_normalized.txt"

    if [[ ! -s "$SFS_OUT/${POP}_sfs_normalized.txt" ]]; then
        echo "WARNING: no data in normalized SFS for $POP"
        continue
    fi

    awk -v pop="$POP" -v n="$N_IND" 'BEGIN{sum=0; singletons=0; doubletons=0}
    NR==1 {singletons=$1}
    NR==2 {doubletons=$1}
    {sum+=$1}
    END {
        if (sum > 0) {
            print "Population:", pop
            print "  Sample size (diploid, folded):", n, "individuals"
            print "  Total segregating sites (sum of folded classes):", sum
            print "  Singletons/lowest-freq class:", singletons
            print "  Doubletons/next class:", doubletons
            print "  Singleton proportion:", singletons/sum
            print "  Doubleton proportion:", doubletons/sum
            print ""
        }
    }' "$SFS_OUT/${POP}_sfs_normalized.txt" >> "$SUMMARY"
done

echo ""
echo "=================================================================="
echo "SFS Analysis Complete"
echo "=================================================================="
cat "$SUMMARY"
echo ""
echo "Individual SFS files in: $SFS_OUT/"
echo "Summary: $SUMMARY"
