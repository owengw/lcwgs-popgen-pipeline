#!/bin/bash
#SBATCH --job-name=p06e_fst_prep
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/p06e_fst_prep_%j.log

# =============================================================================
# p06e_fst_prep.sh  — enumerates all FST comparisons into a job list, then
# submits p06e_fst_array.sh as an array (one task per pair). Replaces the
# serial pairwise loop in the old p06e_fst.sh.
#
# Scales to ANY number of populations: with P pops it builds P(P-1)/2 pairwise
# jobs and runs them 50-at-a-time instead of sequentially.
#
# Usage: sbatch scripts/p06e_fst_prep.sh          (uses OUT_DIR default)
#        OUT_DIR=... sbatch scripts/p06e_fst_prep.sh
# =============================================================================
source ~/.bash_profile
set -eu

OUT_DIR=${OUT_DIR:-$PWD/physalia/angsd_results}
TEMPORAL_FST_THRESHOLD=${TEMPORAL_FST_THRESHOLD:-0.015}
JOBLIST="$OUT_DIR/fst_jobs.tsv"          # one comparison per line: LIST1<TAB>LIST2<TAB>LABEL
mkdir -p logs "$OUT_DIR"
: > "$JOBLIST"

echo "Enumerating FST comparisons -> $JOBLIST"

# --- temporal-significance gate (same logic as before) ---------------------
SKIP_YEAR=false
GLOBAL_FST_FILE="$OUT_DIR/fst_2023_2024.fst.windowed.tsv"
if [[ -f "$GLOBAL_FST_FILE" ]]; then
    G=$(awk 'NR>1 && $5!="nan" && $5!="" {s+=$5;n++} END{if(n>0)printf "%.6f",s/n; else printf "NA"}' "$GLOBAL_FST_FILE")
    if [[ "$G" != "NA" ]]; then
        SIG=$(awk -v f="$G" -v t="$TEMPORAL_FST_THRESHOLD" 'BEGIN{print (f>t)?"YES":"NO"}')
        echo "Global temporal FST=$G  threshold=$TEMPORAL_FST_THRESHOLD  significant=$SIG"
        [[ "$SIG" == "NO" ]] && SKIP_YEAR=true
    fi
else
    echo "NOTE: no global temporal FST yet; including popyear comparisons by default."
fi

# helper: append a comparison if both SAF-able lists exist and aren't already done
add_job() {
    local L1=$1 L2=$2 LABEL=$3
    [[ -f "$L1" && -f "$L2" ]] || return 0
    # n<2 guard: FST/theta are undefined for single-sample populations
    local N1 N2; N1=$(wc -l < "$L1"); N2=$(wc -l < "$L2")
    if (( N1 < 2 || N2 < 2 )); then
        echo "  SKIP $(basename "$L1" .list) vs $(basename "$L2" .list): n<2 (N1=$N1 N2=$N2)"
        return 0
    fi
    local P1 P2; P1=$(basename "$L1" .list); P2=$(basename "$L2" .list)
    [[ -f "$OUT_DIR/fst_${P1}_${P2}.fst.windowed.tsv" ]] && return 0   # skip completed
    printf "%s\t%s\t%s\n" "$L1" "$L2" "$LABEL" >> "$JOBLIST"
}

# --- 1. population-level pairs (auto-discovered, scales to any P) -----------
mapfile -t POP_LISTS < <(ls "$OUT_DIR"/pop_map_lists/pop_map_*.list 2>/dev/null | grep -v "pop_map_ALL.list" | sort)
for ((i=0;i<${#POP_LISTS[@]};i++)); do
  for ((j=i+1;j<${#POP_LISTS[@]};j++)); do
    add_job "${POP_LISTS[i]}" "${POP_LISTS[j]}" "population"
  done
done

# --- 2. within-population temporal -----------------------------------------
for L in "$OUT_DIR"/site_year_map_lists/*2023.list; do
    [[ -e "$L" ]] || continue
    POP=$(basename "$L" 2023.list)
    add_job "$L" "$OUT_DIR/site_year_map_lists/${POP}2024.list" "within-pop-temporal"
done

# --- 3. global temporal -----------------------------------------------------
add_job "$OUT_DIR/year_map_lists/2023.list" "$OUT_DIR/year_map_lists/2024.list" "global-temporal"

# --- 4. between-popyear (only if temporal significant) ---------------------
if [[ "$SKIP_YEAR" == "false" ]]; then
    mapfile -t YL < <(find "$OUT_DIR/site_year_map_lists" -name "*202[34].list" 2>/dev/null | sort)
    for ((i=0;i<${#YL[@]};i++)); do
      for ((j=i+1;j<${#YL[@]};j++)); do
        P1=$(basename "${YL[i]}" .list | sed 's/[0-9]*$//')
        P2=$(basename "${YL[j]}" .list | sed 's/[0-9]*$//')
        [[ "$P1" != "$P2" ]] && add_job "${YL[i]}" "${YL[j]}" "between-popyear"
      done
    done
fi

N=$(wc -l < "$JOBLIST")
echo "Enumerated $N comparisons."
[[ "$N" -eq 0 ]] && { echo "Nothing to do (all done or no lists found)."; exit 0; }

# --- submit the array (50 concurrent) --------------------------------------
AID=$(sbatch --parsable --array=1-${N}%50 \
      --export=ALL,OUT_DIR="$OUT_DIR",JOBLIST="$JOBLIST" \
      scripts/p06e_fst_array.sh)
echo "Submitted FST array $AID ($N tasks, 50 concurrent)."
echo "Collate after it finishes with p06f_summary.sh (or --dependency=afterok:$AID)."
