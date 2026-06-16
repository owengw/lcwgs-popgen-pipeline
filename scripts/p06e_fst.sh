#!/bin/bash
#SBATCH --job-name=p06e_fst
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/p06e_fst_%j.log

# FST calculation script that:
# 1. Checks global temporal FST (2023 vs 2024)
# 2. Only does pop×year comparisons if global FST is significant
# 3. Always does population-level comparisons

source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

set -euo pipefail

OUT_DIR=${OUT_DIR:-$PWD/physalia/angsd_results}
TEMPORAL_FST_THRESHOLD=${TEMPORAL_FST_THRESHOLD:-0.015}
NTHREADS=${SLURM_CPUS_PER_TASK:-8}

echo "=================================================================="
echo "FST Analysis - Temporal Assessment"
echo "=================================================================="

# Check global temporal FST
GLOBAL_FST_FILE="$OUT_DIR/fst_2023_2024.fst.windowed.tsv"

if [[ -f "$GLOBAL_FST_FILE" ]]; then
    GLOBAL_TEMPORAL_FST=$(awk 'NR>1 && $5!="nan" && $5!="" {
        sum += $5; n++
    } END {
        if(n>0) printf "%.6f", sum/n; else printf "NA"
    }' "$GLOBAL_FST_FILE")
    
    if [[ "$GLOBAL_TEMPORAL_FST" != "NA" ]]; then
        IS_GLOBAL_SIG=$(awk -v fst="$GLOBAL_TEMPORAL_FST" -v thresh="$TEMPORAL_FST_THRESHOLD" \
                        'BEGIN {print (fst > thresh) ? "YES" : "NO"}')
        
        echo "Global temporal FST (2023 vs 2024): $GLOBAL_TEMPORAL_FST"
        echo "Threshold: $TEMPORAL_FST_THRESHOLD"
        echo "Significant: $IS_GLOBAL_SIG"
    else
        echo "ERROR: Could not calculate global temporal FST"
        IS_GLOBAL_SIG="UNKNOWN"
    fi
else
    echo "WARNING: Global FST file not found: $GLOBAL_FST_FILE"
    echo "Run FST between 2023 and 2024 first"
    IS_GLOBAL_SIG="UNKNOWN"
fi

# Report per-population temporal FST
echo ""
echo "Per-population temporal FST:"
POPS_WITH_YEARS=("CAI" "NS" "SB" "SI")

for POP in "${POPS_WITH_YEARS[@]}"; do
    FST_FILE="$OUT_DIR/fst_${POP}2023_${POP}2024.fst.windowed.tsv"
    
    if [[ -f "$FST_FILE" ]]; then
        POP_FST=$(awk 'NR>1 && $5!="nan" && $5!="" {
            sum += $5; n++
        } END {
            if(n>0) printf "%.6f", sum/n; else printf "NA"
        }' "$FST_FILE")
        
        if [[ "$POP_FST" != "NA" ]]; then
            echo "  $POP: $POP_FST"
        fi
    fi
done

echo ""

# Decision
if [[ "$IS_GLOBAL_SIG" == "YES" ]]; then
    echo "??  SIGNIFICANT temporal variation (FST=$GLOBAL_TEMPORAL_FST > $TEMPORAL_FST_THRESHOLD)"
    echo "    Will calculate pop×year comparisons"
    SKIP_YEAR_COMPARISONS=false
elif [[ "$IS_GLOBAL_SIG" == "NO" ]]; then
    echo "? No significant temporal variation (FST=$GLOBAL_TEMPORAL_FST = $TEMPORAL_FST_THRESHOLD)"
    echo "  Skipping pop×year comparisons"
    SKIP_YEAR_COMPARISONS=true
else
    echo "??  Cannot determine - defaulting to full analysis"
    SKIP_YEAR_COMPARISONS=false
fi

echo ""
echo "=================================================================="
echo "Running FST Calculations"
echo "=================================================================="
echo ""

# Function to run a single FST comparison
run_fst() {
    local LIST1=$1
    local LIST2=$2
    local LABEL=$3
    
    local POP1=$(basename "$LIST1" .list)
    local POP2=$(basename "$LIST2" .list)
    local OUTNAME="${POP1}_${POP2}"
    
    # Skip if already done
    if [[ -f "$OUT_DIR/fst_${OUTNAME}.fst.windowed.tsv" ]]; then
        echo "  ? $OUTNAME (already exists)"
        return 0
    fi
    
    echo "  Running: $OUTNAME ($LABEL)"
    
    # Get SAF files
    local SAF1="$OUT_DIR/${POP1}.saf.idx"
    local SAF2="$OUT_DIR/${POP2}.saf.idx"
    
    if [[ ! -f "$SAF1" ]] || [[ ! -f "$SAF2" ]]; then
        echo "    WARNING: SAF files not found, skipping"
        return 1
    fi
    
    # Calculate 2D SFS
    realSFS "$SAF1" "$SAF2" -P $NTHREADS > "$OUT_DIR/fst_${OUTNAME}.2dsfs" 2>/dev/null
    
    # Calculate FST index
    realSFS fst index "$SAF1" "$SAF2" \
        -sfs "$OUT_DIR/fst_${OUTNAME}.2dsfs" \
        -fstout "$OUT_DIR/fst_${OUTNAME}" \
        -P $NTHREADS 2>/dev/null
    
    # Get global FST
    realSFS fst stats "$OUT_DIR/fst_${OUTNAME}.fst.idx" \
        > "$OUT_DIR/fst_${OUTNAME}.fst.global.txt" 2>/dev/null
    
    # Get windowed FST (50kb windows, 10kb step)
    realSFS fst stats2 "$OUT_DIR/fst_${OUTNAME}.fst.idx" \
        -win 50000 -step 10000 \
        > "$OUT_DIR/fst_${OUTNAME}.fst.windowed.tsv" 2>/dev/null
    
    echo "    ? Complete"
}

# Counter
TOTAL=0
COMPLETED=0

# 1. Population-level comparisons (ALWAYS)
echo "Population-level comparisons:"
POPULATIONS=("pop_map_CAI" "pop_map_COI" "pop_map_NS" "pop_map_SB" "pop_map_SI" "pop_map_SP")

for ((i=0; i<${#POPULATIONS[@]}; i++)); do
    for ((j=i+1; j<${#POPULATIONS[@]}; j++)); do
        LIST1="$OUT_DIR/pop_map_lists/${POPULATIONS[i]}.list"
        LIST2="$OUT_DIR/pop_map_lists/${POPULATIONS[j]}.list"
        
        if [[ -f "$LIST1" && -f "$LIST2" ]]; then
            ((TOTAL++))
            if run_fst "$LIST1" "$LIST2" "population"; then
                ((COMPLETED++))
            fi
        fi
    done
done

# 2. Within-population temporal (ALWAYS - for reference)
echo ""
echo "Within-population temporal comparisons:"
for POP in "${POPS_WITH_YEARS[@]}"; do
    LIST1="$OUT_DIR/site_year_map_lists/${POP}2023.list"
    LIST2="$OUT_DIR/site_year_map_lists/${POP}2024.list"
    
    [[ ! -f "$LIST1" ]] && LIST1="$OUT_DIR/pop_map_lists/${POP}2023.list"
    [[ ! -f "$LIST2" ]] && LIST2="$OUT_DIR/pop_map_lists/${POP}2024.list"
    
    if [[ -f "$LIST1" && -f "$LIST2" ]]; then
        ((TOTAL++))
        if run_fst "$LIST1" "$LIST2" "within-pop temporal"; then
            ((COMPLETED++))
        fi
    fi
done

# 3. Global temporal (ALWAYS - for reference)
echo ""
echo "Global temporal comparison:"
YEAR_2023="$OUT_DIR/year_map_lists/2023.list"
YEAR_2024="$OUT_DIR/year_map_lists/2024.list"
[[ ! -f "$YEAR_2023" ]] && YEAR_2023="$OUT_DIR/pop_map_lists/2023.list"
[[ ! -f "$YEAR_2024" ]] && YEAR_2024="$OUT_DIR/pop_map_lists/2024.list"

if [[ -f "$YEAR_2023" && -f "$YEAR_2024" ]]; then
    ((TOTAL++))
    if run_fst "$YEAR_2023" "$YEAR_2024" "global temporal"; then
        ((COMPLETED++))
    fi
fi

# 4. Between-population×year (CONDITIONAL)
if [[ "$SKIP_YEAR_COMPARISONS" == "false" ]]; then
    echo ""
    echo "Between-population×year comparisons:"
    
    # Get year-specific lists
    YEAR_LISTS=($(find $OUT_DIR/site_year_map_lists -name "*202[34].list" 2>/dev/null | sort))
    [[ ${#YEAR_LISTS[@]} -eq 0 ]] && YEAR_LISTS=($(find $OUT_DIR/pop_map_lists -name "*202[34].list" ! -name "2023.list" ! -name "2024.list" 2>/dev/null | sort))
    
    for ((i=0; i<${#YEAR_LISTS[@]}; i++)); do
        for ((j=i+1; j<${#YEAR_LISTS[@]}; j++)); do
            LIST1="${YEAR_LISTS[i]}"
            LIST2="${YEAR_LISTS[j]}"
            
            BASE1=$(basename "$LIST1" .list)
            BASE2=$(basename "$LIST2" .list)
            POP1=$(echo "$BASE1" | sed 's/[0-9]*$//')
            POP2=$(echo "$BASE2" | sed 's/[0-9]*$//')
            
            # Skip if same population
            if [[ "$POP1" != "$POP2" ]]; then
                ((TOTAL++))
                if run_fst "$LIST1" "$LIST2" "between-pop×year"; then
                    ((COMPLETED++))
                fi
            fi
        done
    done
else
    echo ""
    echo "Skipping between-population×year comparisons (temporal variation not significant)"
fi

echo ""
echo "=================================================================="
echo "FST Analysis Complete"
echo "=================================================================="
echo "Total comparisons: $TOTAL"
echo "Completed: $COMPLETED"
echo "Already existed: $((TOTAL - COMPLETED))"
echo ""