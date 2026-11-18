#!/bin/bash
#SBATCH --job-name=p06c_angsd
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=logs/p06c_angsd_%A.log

# EDITS:
#  - set -eu (no pipefail)
#  - CHROMOSOME FILTER: restrict ANGSD to the 13 RagTag chromosomes via -rf,
#    excluding unplaced contig_*_RagTag_pilon junk. Controlled by CHROM_FILE
#    (default physalia/angsd_results/chromosomes.txt). Kept the existing
#    random-sites option; the two combine correctly.
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/owenspopgen || true
set -eu
set -x

ANGSD_BIN=angsd
BAM_LIST=${BAM_LIST:-}
REFERENCE=${REFERENCE:-}
OUT_DIR=${OUT_DIR:-}
ANCESTRAL=${ANCESTRAL:-}
THREADS=${THREADS:-8}
MIN_MAPQ=${MIN_MAPQ:-10}
MIN_Q=${MIN_Q:-13}
FOLD_SFS=${FOLD_SFS:-1}
MIN_IND_RATIO=${MIN_IND_RATIO:-0.20}
USE_RANDOM_SITES=${USE_RANDOM_SITES:-0}
# NEW: chromosome restriction (set CHROM_FILE="" to disable)
CHROM_FILE=${CHROM_FILE:-physalia/angsd_results/chromosomes.txt}

[[ -f "$BAM_LIST" ]] || { echo "ERROR: BAM list not found: $BAM_LIST"; exit 1; }
[[ -f "$REFERENCE.fai" ]] || { echo "ERROR: Reference FAI not found: $REFERENCE.fai"; exit 1; }

POP_NAME=$(basename "$BAM_LIST" .list)
OUT_PREFIX="$OUT_DIR/$POP_NAME"
mkdir -p "$OUT_DIR/logs"

N_SAMPLES=$(wc -l < "$BAM_LIST")
MIN_IND=$(echo "$N_SAMPLES * $MIN_IND_RATIO" | bc | awk '{print int($1+0.5)}')
[[ $MIN_IND -lt 1 ]] && MIN_IND=1
[[ $MIN_IND -gt $N_SAMPLES ]] && MIN_IND=$N_SAMPLES

echo "=================================================================="
echo "ANGSD SAF: $POP_NAME | samples $N_SAMPLES | minInd $MIN_IND"
[[ -n "$CHROM_FILE" && -f "$CHROM_FILE" ]] && echo "Chromosome filter: $CHROM_FILE ($(wc -l < "$CHROM_FILE") chromosomes)"
echo "=================================================================="

if [[ -n "$ANCESTRAL" && -f "$ANCESTRAL" ]]; then ANC_OPT="-anc $ANCESTRAL"; else ANC_OPT="-anc $REFERENCE"; fi

SITES_OPT=""
RF_OPT=""

# chromosome restriction (default). If random-sites is on, that path sets its
# own -rf (from the sites' contigs) which already lie within chromosomes.
if [[ -n "$CHROM_FILE" && -f "$CHROM_FILE" ]]; then
    RF_OPT="-rf $CHROM_FILE"
fi

if [[ "$USE_RANDOM_SITES" == "1" ]]; then
    if [[ ! -f "$OUT_DIR/random_10M_sites.list" ]]; then
        echo "Generating random 10M sites from chromosomes >100kb..."
        # restrict site generation to chromosomes only
        awk 'NR==FNR{keep[$1]=1; next} ($1 in keep) && $2>100000' \
            "$CHROM_FILE" "$REFERENCE.fai" | \
            awk 'BEGIN{srand(42)}{for(i=1;i<=$2;i+=150) if(rand()<0.067) print $1"\t"i"\tA\tT"}' \
            | head -10000000 > "$OUT_DIR/random_10M_sites.list"
        angsd sites index "$OUT_DIR/random_10M_sites.list"
        cut -f1 "$OUT_DIR/random_10M_sites.list" | sort -u > "$OUT_DIR/random_10M_contigs.txt"
    fi
    SITES_OPT="-sites $OUT_DIR/random_10M_sites.list"
    RF_OPT="-rf $OUT_DIR/random_10M_contigs.txt"
fi

echo "Running ANGSD SAF..."
if ! $ANGSD_BIN -bam "$BAM_LIST" -ref "$REFERENCE" $ANC_OPT \
    -GL 2 -doSaf 1 -doCounts 1 \
    -minMapQ $MIN_MAPQ -minQ $MIN_Q -minInd $MIN_IND \
    -setMinDepthInd 1 -setMaxDepthInd 50 -P $THREADS \
    -nQueueSize 50 -remove_bads 1 -trim 0 \
    $SITES_OPT $RF_OPT \
    -out "$OUT_PREFIX" 2>&1 | tee "$OUT_PREFIX.angsd_saf.log"; then
    echo "ERROR: ANGSD SAF failed"; tail -50 "$OUT_PREFIX.angsd_saf.log"; exit 1
fi

[[ -f "$OUT_PREFIX.saf.idx" ]] || { echo "ERROR: SAF index not created"; exit 1; }
echo "SAF complete: $OUT_PREFIX.saf.idx"
