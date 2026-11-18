#!/bin/bash
#SBATCH --job-name=theta_to_corrected
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/p06_theta_to_corrected_%j.log

# =============================================================================
# p06_theta_to_corrected.sh — ADAPTER, replaces p07b_pop_theta.sh entirely.
#
# p06 already computed all-sites, per-chromosome, genome-wide theta for
# CAI/NS/SB/SI (pop_map_<POP>.pestPG, via p06d_theta_bychrom.sh --merge).
# Re-running p07b's own SAF->SFS->saf2theta pipeline on random-10M sites
# would be both redundant AND a downgrade (sparser data than what we have).
#
# This script just reformats the existing pestPG files into the exact
# theta_corrected/theta_summary_corrected.tsv format that p07c_inbreeding.sh
# and the local R plotting scripts (plot_diversity.R, pop_gen_analysis.R)
# expect, so nothing downstream needs touching.
#
# pestPG columns (thetaStat do_stat header):
#   1:indices 2:Chr 3:WinCenter 4:tW 5:tP 6:tF 7:tH 8:tL 9:Tajima
#   10:fuf 11:fud 12:fayh 13:zeng 14:nSites
# =============================================================================

source ~/.bash_profile
set -eu

OUT_DIR=${OUT_DIR:-physalia/angsd_results}
THETA_DIR="$OUT_DIR/theta_corrected"
POPULATIONS=(CAI NS SB SI)   # DH/COI/SP excluded — matches repo convention (n too small)

mkdir -p "$THETA_DIR"
SUMMARY="$THETA_DIR/theta_summary_corrected.tsv"
echo -e "population\tn_sites\tWatterson_theta\ttheta_pi\ttajima_D\tWatterson_theta_per_site\ttheta_pi_per_site" > "$SUMMARY"

for POP in "${POPULATIONS[@]}"; do
    SRC="$OUT_DIR/pop_map_${POP}.pestPG"
    if [[ ! -s "$SRC" ]]; then
        echo "WARNING: $SRC not found — skipping $POP (has its p06d theta chain finished?)"
        continue
    fi

    STATS=$(awk -F'\t' 'NR>1 && NF>=14 {
        tw_sum += $4; tp_sum += $5; taj_sum += $9; sites_sum += $14; n_windows++
    } END {
        if (n_windows > 0 && sites_sum > 0) {
            mean_tw  = tw_sum  / n_windows
            mean_tp  = tp_sum  / n_windows
            mean_taj = taj_sum / n_windows
            printf "%d\t%.6f\t%.6f\t%.6f\t%.10f\t%.10f",
                   sites_sum, mean_tw, mean_tp, mean_taj,
                   tw_sum/sites_sum, tp_sum/sites_sum
        } else {
            printf "0\tNA\tNA\tNA\tNA\tNA"
        }
    }' "$SRC")

    echo -e "pop_map_${POP}\t$STATS" >> "$SUMMARY"
    echo "$POP: $STATS"
done

echo ""
echo "=================================================================="
echo "theta_summary_corrected.tsv written -> $SUMMARY"
column -t "$SUMMARY"
echo "=================================================================="
