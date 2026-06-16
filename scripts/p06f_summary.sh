#!/bin/bash
#SBATCH --job-name=p06f_summary
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/p06f_summary_%A.log

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -euo pipefail

OUT_DIR=${OUT_DIR:-$PWD/physalia/angsd_results}

echo "=================================================================="
echo "Generating Summary Statistics"
echo "=================================================================="

# -------------------------
# THETA SUMMARY (CORRECTED)
# -------------------------
THETA_SUMMARY="$OUT_DIR/theta_summary.tsv"

echo "Generating theta summary..."
echo -e "population\tn_sites\tWatterson_theta\ttheta_pi\ttajima_D\tWatterson_theta_per_site\ttheta_pi_per_site" > "$THETA_SUMMARY"

for FILE in "$OUT_DIR"/*.thetas.windowed.tsv; do
    [[ -f "$FILE" ]] || continue

    GROUP=$(basename "$FILE" .thetas.windowed.tsv)

    # thetaStat windowed output columns:
    # (indexStart,indexStop),(firstPos_withData,lastPos_withData),WinCenter,tW,tP,tF,tH,tL,Tajima,fuf,fud,fayh,zeng,nSites
    # We want: tW (col 4), tP (col 5), Tajima (col 9), nSites (col 14)
    
    # Calculate TOTAL and MEAN per-site
    STATS=$(awk -F'\t' 'NR>1 && NF>=14 && $1!~/^#/ && $1!~/^Chr/ {
        tw_sum += $4
        tp_sum += $5
        taj_sum += $9
        sites_sum += $14
        n_windows++
    } END {
        if(n_windows > 0 && sites_sum > 0) {
            mean_tw = tw_sum / n_windows
            mean_tp = tp_sum / n_windows
            mean_taj = taj_sum / n_windows
            tw_per_site = tw_sum / sites_sum
            tp_per_site = tp_sum / sites_sum
            printf "%d\t%.6f\t%.6f\t%.6f\t%.10f\t%.10f", sites_sum, mean_tw, mean_tp, mean_taj, tw_per_site, tp_per_site
        } else {
            printf "0\tNA\tNA\tNA\tNA\tNA"
        }
    }' "$FILE")

    echo -e "$GROUP\t$STATS" >> "$THETA_SUMMARY"
done

echo "? Theta summary created: $THETA_SUMMARY"

# -------------------------
# FST SUMMARY (CORRECTED)
# -------------------------
FST_SUMMARY="$OUT_DIR/fst_summary.tsv"

echo "Generating FST summary..."
echo -e "Pop1\tPop2\tFST_mean\tFST_weighted\tNum_Windows" > "$FST_SUMMARY"

for FILE in "$OUT_DIR"/fst_*.fst.windowed.tsv; do
    [[ -f "$FILE" ]] || continue
    
    BASE=$(basename "$FILE" .fst.windowed.tsv)
    # Remove 'fst_' prefix
    POPS=$(echo "$BASE" | sed 's/^fst_//')
    
    # Extract POP1 and POP2 using regex patterns to handle underscores in names
    if [[ "$POPS" =~ ^(.+)_(pop_map_.+)$ ]]; then
        # Handles: pop_map_CAI_pop_map_NS
        POP1="${BASH_REMATCH[1]}"
        POP2="${BASH_REMATCH[2]}"
    elif [[ "$POPS" =~ ^(pop_map_[A-Z]+)_([A-Z]+[0-9]{4})$ ]]; then
        # Handles: pop_map_CAI_NS2023
        POP1="${BASH_REMATCH[1]}"
        POP2="${BASH_REMATCH[2]}"
    elif [[ "$POPS" =~ ^([A-Z]+[0-9]{4})_(pop_map_[A-Z]+)$ ]]; then
        # Handles: CAI2023_pop_map_NS
        POP1="${BASH_REMATCH[1]}"
        POP2="${BASH_REMATCH[2]}"
    elif [[ "$POPS" =~ ^([A-Z]+[0-9]{4})_([A-Z]+[0-9]{4})$ ]]; then
        # Handles: CAI2023_NS2024
        POP1="${BASH_REMATCH[1]}"
        POP2="${BASH_REMATCH[2]}"
    elif [[ "$POPS" =~ ^([0-9]{4})_([0-9]{4})$ ]]; then
        # Handles: 2023_2024
        POP1="${BASH_REMATCH[1]}"
        POP2="${BASH_REMATCH[2]}"
    else
        # Fallback
        POP1=$(echo "$POPS" | cut -d'_' -f1)
        POP2=$(echo "$POPS" | sed 's/^[^_]*_//')
    fi
    
    # Calculate mean and weighted FST
    # thetaStat FST output columns: region chr midPos Nsites fst
    STATS=$(awk 'NR>1 && $5!="nan" && $5!="NA" && $5!="" {
        sum_fst += $5
        sum_weighted += $5 * $4
        sum_sites += $4
        n++
    } END {
        if(n>0 && sum_sites>0) {
            mean_fst = sum_fst / n
            weighted_fst = sum_weighted / sum_sites
            printf "%.6f\t%.6f\t%d", mean_fst, weighted_fst, n
        } else {
            printf "NA\tNA\t0"
        }
    }' "$FILE")
    
    echo -e "$POP1\t$POP2\t$STATS" >> "$FST_SUMMARY"
done

echo "? FST summary created: $FST_SUMMARY"

echo ""
echo "=================================================================="
echo "Summary Statistics Complete"
echo "=================================================================="
echo "Files created:"
echo "  - $THETA_SUMMARY"
echo "  - $FST_SUMMARY"
echo "=================================================================="