#!/bin/bash
#SBATCH --job-name=theta_downsampled
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/p07d_theta_downsampled_%j.log

# =============================================================================
# p07d_theta_downsampled.sh — RagTag rewrite. Equal-n theta to remove sample
# size bias from Tajima's D comparisons across CAI/NS/SB/SI.
#
# CHANGE FROM REPO DEFAULT: DOWNSAMPLE_N=17 assumed SB=17. Your current SB
# is 16 (post SI35/replicate exclusion), so DOWNSAMPLE_N=16 — the smallest
# of the four main populations. SB is used as-is (all 16); CAI/NS/SI are
# randomly subsampled to 16, seed=42.
#
# This script does NOT reimplement SAF/theta itself — it builds the
# subsampled BAM lists, then DISPATCHES to the already-tested per-chromosome
# infrastructure (p06c_saf_bychrom.sh + p06d_theta_bychrom.sh) with a
# "_ds16" population suffix, so it inherits the OOM/timeout fixes already
# hard-won this session rather than repeating them. Submits jobs and exits;
# does not block waiting for them (SLURM chain via --dependency).
# =============================================================================

source ~/.bash_profile
set -eu

OUT_DIR=${OUT_DIR:-physalia/angsd_results}
REFERENCE=${REFERENCE:-/mnt/parscratch/users/bi4og/ncbi_submission/plongirostris.fasta}
POP_LIST_DIR="$OUT_DIR/pop_map_lists"
DS_DIR="$OUT_DIR/theta_downsampled"
DOWNSAMPLE_N=16
SEED=42
NCHR=$(wc -l < "$OUT_DIR/chromosomes.txt")

POPULATIONS=(CAI NS SB SI)

mkdir -p "$DS_DIR/subsampled_lists" logs

echo "=================================================================="
echo "Building downsampled BAM lists (n=$DOWNSAMPLE_N, seed=$SEED)"
echo "=================================================================="

for POP in "${POPULATIONS[@]}"; do
    FULL_LIST="$POP_LIST_DIR/pop_map_${POP}.list"
    SUB_LIST="$DS_DIR/subsampled_lists/${POP}_ds${DOWNSAMPLE_N}.list"

    [[ -f "$FULL_LIST" ]] || { echo "ERROR: $FULL_LIST not found"; exit 1; }
    N_FULL=$(wc -l < "$FULL_LIST")

    if [[ -s "$SUB_LIST" ]]; then
        echo "  $POP: subsampled list already exists ($(wc -l < "$SUB_LIST") lines), skipping"
        continue
    fi

    if [[ $N_FULL -lt $DOWNSAMPLE_N ]]; then
        echo "ERROR: $POP has only $N_FULL individuals, need $DOWNSAMPLE_N"; exit 1
    elif [[ $N_FULL -eq $DOWNSAMPLE_N ]]; then
        cp "$FULL_LIST" "$SUB_LIST"
        echo "  $POP: already at n=$DOWNSAMPLE_N, using full list"
    else
        shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:"seed${SEED}" -nosalt </dev/zero 2>/dev/null) \
            -n $DOWNSAMPLE_N "$FULL_LIST" > "$SUB_LIST"
        echo "  $POP: subsampled $N_FULL -> $DOWNSAMPLE_N"
    fi
done

echo ""
echo "=================================================================="
echo "Dispatching SAF -> SFS -> theta chain per population (ds${DOWNSAMPLE_N})"
echo "=================================================================="

for POP in "${POPULATIONS[@]}"; do
    DSPOP="pop_map_${POP}_ds${DOWNSAMPLE_N}"
    SUB_LIST="$DS_DIR/subsampled_lists/${POP}_ds${DOWNSAMPLE_N}.list"

    if [[ -s "$OUT_DIR/${DSPOP}.pestPG" ]]; then
        echo "  $DSPOP: theta already complete, skipping"
        continue
    fi

    # 1. per-chromosome SAF array
    AID=$(sbatch --parsable --array=1-$NCHR%13 --job-name="saf_${DSPOP}" \
        --export=ALL,POP="$DSPOP",BAM_LIST="$PWD/$SUB_LIST",REFERENCE="$REFERENCE",OUT_DIR="$PWD/$OUT_DIR" \
        scripts/p06c_saf_bychrom.sh)
    # 2. merge SAF parts
    MID=$(sbatch --parsable --dependency=afterok:$AID --mem=16G --time=04:00:00 --job-name="safcat_${DSPOP}" \
        --wrap="bash scripts/p06c_saf_bychrom.sh --merge POP=$DSPOP OUT_DIR=$PWD/$OUT_DIR")
    # 3. per-chromosome SFS array (on the merged SAF's parts — reuse saf_parts/ directly)
    SID=$(sbatch --parsable --dependency=afterok:$MID --array=1-$NCHR%13 --job-name="sfs_${DSPOP}" \
        --export=ALL,POP="$DSPOP",OUT_DIR="$PWD/$OUT_DIR" \
        scripts/p06d_theta_bychrom.sh --sfs-array)
    # 4. sum SFS
    SUMID=$(sbatch --parsable --dependency=afterok:$SID --job-name="sfssum_${DSPOP}" \
        --wrap="bash scripts/p06d_theta_bychrom.sh --sfs-sum POP=$DSPOP OUT_DIR=$PWD/$OUT_DIR")
    # 5. per-chromosome theta array
    TID=$(sbatch --parsable --dependency=afterok:$SUMID --array=1-$NCHR%13 --time=16:00:00 --job-name="theta_${DSPOP}" \
        --export=ALL,POP="$DSPOP",OUT_DIR="$PWD/$OUT_DIR" \
        scripts/p06d_theta_bychrom.sh)
    # 6. merge theta
    sbatch --dependency=afterok:$TID --job-name="thetamerge_${DSPOP}" \
        --wrap="bash scripts/p06d_theta_bychrom.sh --merge POP=$DSPOP OUT_DIR=$PWD/$OUT_DIR"

    echo "  $DSPOP: chain submitted (SAF=$AID merge=$MID sfs=$SID sum=$SUMID theta=$TID)"
done

echo ""
echo "Once all chains complete, collate with:"
echo "  bash scripts/p07d_collate_downsampled.sh"
