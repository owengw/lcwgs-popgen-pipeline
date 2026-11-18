#!/bin/bash
#SBATCH --job-name=p06e_fst
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/p06e_fst_%A_%a.log

# =============================================================================
# p06e_fst_array.sh — runs ONE FST comparison (the SLURM_ARRAY_TASK_ID-th line
# of $JOBLIST). Submitted by p06e_fst_prep.sh. set -eu (no pipefail).
# =============================================================================
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen
set -eu

OUT_DIR=${OUT_DIR:-$PWD/physalia/angsd_results}
JOBLIST=${JOBLIST:-$OUT_DIR/fst_jobs.tsv}
NTHREADS=${SLURM_CPUS_PER_TASK:-8}

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
[[ -n "$LINE" ]] || { echo "ERROR: no job on line $SLURM_ARRAY_TASK_ID"; exit 1; }
LIST1=$(echo "$LINE" | cut -f1); LIST2=$(echo "$LINE" | cut -f2); LABEL=$(echo "$LINE" | cut -f3)

POP1=$(basename "$LIST1" .list); POP2=$(basename "$LIST2" .list)
OUTNAME="${POP1}_${POP2}"
SAF1="$OUT_DIR/${POP1}.saf.idx"; SAF2="$OUT_DIR/${POP2}.saf.idx"

echo "FST: $OUTNAME ($LABEL)"
[[ -f "$SAF1" && -f "$SAF2" ]] || { echo "ERROR: SAF missing ($SAF1 / $SAF2) - run p06c first"; exit 1; }
[[ -f "$OUT_DIR/fst_${OUTNAME}.fst.windowed.tsv" ]] && { echo "already done, skipping"; exit 0; }

# realSFS: stdout MUST go straight to file (no tee) or the SFS is corrupted.
realSFS "$SAF1" "$SAF2" -P $NTHREADS > "$OUT_DIR/fst_${OUTNAME}.2dsfs" 2>"$OUT_DIR/fst_${OUTNAME}.2dsfs.log"
realSFS fst index "$SAF1" "$SAF2" -sfs "$OUT_DIR/fst_${OUTNAME}.2dsfs" \
    -fstout "$OUT_DIR/fst_${OUTNAME}" -P $NTHREADS 2>"$OUT_DIR/fst_${OUTNAME}.index.log"
realSFS fst stats "$OUT_DIR/fst_${OUTNAME}.fst.idx" > "$OUT_DIR/fst_${OUTNAME}.fst.global.txt" 2>/dev/null
realSFS fst stats2 "$OUT_DIR/fst_${OUTNAME}.fst.idx" -win 50000 -step 25000 \
    > "$OUT_DIR/fst_${OUTNAME}.fst.windowed.tsv" 2>/dev/null

[[ -s "$OUT_DIR/fst_${OUTNAME}.fst.windowed.tsv" ]] && echo "Complete: $OUTNAME" || { echo "ERROR: no windowed output"; exit 1; }
