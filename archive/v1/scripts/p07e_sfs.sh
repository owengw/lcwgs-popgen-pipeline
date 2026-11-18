#!/bin/bash
#SBATCH --job-name=sfs_analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/p07e_sfs_%j.log

# Load environment
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen

# Site Frequency Spectrum analysis and plotting

set -euo pipefail

# Required variables
: ${OUT_DIR:?}

echo "=================================================================="
echo "Site Frequency Spectrum Analysis"
echo "=================================================================="
echo "Output directory: $OUT_DIR"
echo "=================================================================="

mkdir -p ${OUT_DIR}/logs
mkdir -p ${OUT_DIR}/sfs

# Initialize summary file
echo "SFS Summary Statistics" > ${OUT_DIR}/sfs/sfs_summary.txt
echo "======================" >> ${OUT_DIR}/sfs/sfs_summary.txt
echo "" >> ${OUT_DIR}/sfs/sfs_summary.txt

# Find all SFS files
echo "Processing SFS files for all populations..."

for POP_NAME in pop_map_CAI pop_map_COI pop_map_DH pop_map_NS pop_map_SB pop_map_SI pop_map_SP \
                CAI2023 CAI2024 COI2023 DH2023 NS2023 NS2024 SB2023 SB2024 SI2023 SI2024 SP2023 \
                2023 2024 pop_map_ALL; do
    
    SFS_FILE="${OUT_DIR}/${POP_NAME}.sfs"
    
    if [[ ! -f "${SFS_FILE}" ]]; then
        echo "WARNING: SFS file not found for ${POP_NAME}"
        continue
    fi
    
    echo "Processing ${POP_NAME}..."
    
    # Get sample size from the corresponding BAM list
    N_IND=""
    
    # Try different list directories
    for LIST_DIR in pop_map_lists site_year_map_lists year_map_lists; do
        if [[ -f "${OUT_DIR}/${LIST_DIR}/${POP_NAME}.list" ]]; then
            N_IND=$(wc -l < ${OUT_DIR}/${LIST_DIR}/${POP_NAME}.list)
            break
        fi
    done
    
    # If still not found, try to infer from SFS file size
    if [[ -z "$N_IND" ]]; then
        echo "WARNING: Cannot find BAM list for ${POP_NAME}, attempting to infer from SFS"
        # Count number of values in SFS file (should be 2*N + 1)
        N_VALUES=$(cat ${SFS_FILE} | wc -w)
        N_IND=$(( (N_VALUES - 1) / 2 ))
        
        if [[ $N_IND -le 0 ]]; then
            echo "ERROR: Cannot determine sample size for ${POP_NAME}"
            continue
        fi
        echo "  Inferred sample size: ${N_IND}"
    fi
    
    echo "  Sample size: ${N_IND} individuals (${N_IND}*2 = $((N_IND*2)) chromosomes)"
    
    # Create normalized SFS (exclude first and last bins - invariant sites)
    # SFS format: space-separated values for allele frequencies from 0 to 2N
    cat ${SFS_FILE} | tr ' ' '\n' | awk 'NF>0 && NR>1 && NR<='$((N_IND*2))' {print}' > ${OUT_DIR}/sfs/${POP_NAME}_sfs_normalized.txt
    
    # Check if normalized file has content
    if [[ ! -s ${OUT_DIR}/sfs/${POP_NAME}_sfs_normalized.txt ]]; then
        echo "WARNING: No data in normalized SFS for ${POP_NAME}"
        continue
    fi
    
    # Calculate SFS statistics
    awk -v pop="${POP_NAME}" -v n=$((N_IND*2)) 'BEGIN {
        sum=0
        singletons=0
        doubletons=0
    }
    NR==1 {singletons=$1}
    NR==2 {doubletons=$1}
    {sum+=$1}
    END {
        if (sum > 0) {
            print "Population:", pop
            print "  Sample size:", n/2, "individuals"
            print "  Total segregating sites:", sum
            print "  Singletons:", singletons
            print "  Doubletons:", doubletons
            print "  Singleton proportion:", singletons/sum
            print "  Doubleton proportion:", doubletons/sum
            print ""
        }
    }' ${OUT_DIR}/sfs/${POP_NAME}_sfs_normalized.txt >> ${OUT_DIR}/sfs/sfs_summary.txt
    
done

echo ""
echo "=================================================================="
echo "SFS Analysis Complete"
echo "=================================================================="
echo ""
echo "SFS summary:"
cat ${OUT_DIR}/sfs/sfs_summary.txt
echo ""
echo "Individual SFS files in: ${OUT_DIR}/sfs/"
echo "Summary: ${OUT_DIR}/sfs/sfs_summary.txt"
echo "=================================================================="