#!/bin/bash
#SBATCH --job-name=p00_preflight
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/p00_preflight_%j.log

# =============================================================================
# p00_preflight.sh — run ONCE before p06. Makes the reference and sample lists
# safe for the rest of the pipeline, and fails loudly if anything is wrong.
#
# Does:
#   1. Normalises the reference: a REAL .fasta (not a symlink, not .fsa),
#      because Picard rejects .fsa and dangling/symlinked .fai break ANGSD.
#   2. Builds .fai, .dict, and (optionally) the bowtie2 index if missing.
#   3. Writes chromosomes.txt (placed chromosomes only) for the -rf filter.
#   4. Validates the BAM list: files exist, pass quickcheck, and their @SQ
#      names match the reference (the #1 cause of "ANGSD retained 0 sites").
#   5. Clears stale, regenerable files that silently poison a rerun.
#
# Usage:
#   sbatch scripts/p00_preflight.sh -r <reference> -b <bam_list> -o <out_dir> [-c <chr_regex>] [-B]
#     -r  reference FASTA (may be .fsa / a symlink; will be normalised)
#     -b  BAM list (one path per line)
#     -o  ANGSD output dir (default: physalia/angsd_results)
#     -c  regex matching placed chromosome names in the FASTA headers
#         (default: '^chr[0-9]+_RagTag_pilon$')
#     -B  also build the bowtie2 index if absent (slow; only needed pre-alignment)
# =============================================================================

_SAVED_ARGS=("$@")
source ~/.bash_profile
conda activate /mnt/parscratch/users/bi4og/conda_envs/align39   # samtools + picard, glibc-safe
set -- "${_SAVED_ARGS[@]}"
set -eu

REFERENCE=""; BAM_LIST=""; OUT_DIR="physalia/angsd_results"
CHR_REGEX='^chr[0-9]+_RagTag_pilon$'; BUILD_BT2=0

usage(){ sed -n '10,26p' "$0"; exit 1; }
while getopts "r:b:o:c:Bh" opt; do
  case $opt in
    r) REFERENCE="$OPTARG" ;;
    b) BAM_LIST="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    c) CHR_REGEX="$OPTARG" ;;
    B) BUILD_BT2=1 ;;
    *) usage ;;
  esac
done
[[ -n "$REFERENCE" && -n "$BAM_LIST" ]] || { echo "ERROR: -r and -b are required"; usage; }
mkdir -p logs "$OUT_DIR"

fail(){ echo "PREFLIGHT FAILED: $*"; exit 1; }
ok(){   echo "  OK: $*"; }

echo "=================================================================="
echo "1. Normalising reference"
echo "=================================================================="
# force a .fasta extension (Picard/htsjdk reject .fsa)
case "$REFERENCE" in
    *.fasta|*.fa|*.fna) ;;
    *) NEW="${REFERENCE%.*}.fasta"
       echo "  $REFERENCE has a non-standard extension -> using $NEW"
       [[ -e "$NEW" ]] || cp "$REFERENCE" "$NEW"
       REFERENCE="$NEW" ;;
esac
# resolve symlinks into a real file (tools that follow the link then see .fsa)
if [[ -L "$REFERENCE" ]]; then
    TARGET=$(readlink -f "$REFERENCE")
    echo "  $REFERENCE is a symlink -> materialising a real copy of $TARGET"
    cp --remove-destination "$TARGET" "$REFERENCE"
fi
[[ -s "$REFERENCE" ]] || fail "reference missing/empty: $REFERENCE"
ok "reference is a real file: $REFERENCE ($(du -h "$REFERENCE" | cut -f1))"

echo "=================================================================="
echo "2. Indices"
echo "=================================================================="
# .fai — remove dangling/symlinked ones and rebuild
if [[ -L "$REFERENCE.fai" || ! -s "$REFERENCE.fai" ]]; then
    rm -f "$REFERENCE.fai"; samtools faidx "$REFERENCE"
fi
[[ -s "$REFERENCE.fai" && ! -L "$REFERENCE.fai" ]] || fail "could not create real $REFERENCE.fai"
ok ".fai present ($(wc -l < "$REFERENCE.fai") sequences)"

DICT="${REFERENCE%.*}.dict"
if [[ ! -s "$DICT" ]]; then
    rm -f "$DICT"                      # picard refuses to overwrite
    picard CreateSequenceDictionary R="$REFERENCE" O="$DICT"
fi
ok ".dict present"

BT2_BASE="${REFERENCE%.*}"
if [[ "$BUILD_BT2" -eq 1 ]]; then
    if ls "${BT2_BASE}".*.bt2 >/dev/null 2>&1; then
        ok "bowtie2 index already present"
    else
        echo "  building bowtie2 index (slow)..."
        bowtie2-build --threads "${SLURM_CPUS_PER_TASK:-8}" "$REFERENCE" "$BT2_BASE"
    fi
    N_BT2=$(ls "${BT2_BASE}".*.bt2 2>/dev/null | wc -l)
    [[ "$N_BT2" -ge 6 ]] || fail "bowtie2 index incomplete ($N_BT2 files, expected 6)"
    ok "bowtie2 index complete"
fi

echo "=================================================================="
echo "3. Chromosome filter"
echo "=================================================================="
CHROM_FILE="$OUT_DIR/chromosomes.txt"
cut -f1 "$REFERENCE.fai" | grep -E "$CHR_REGEX" > "$CHROM_FILE" || true
[[ -s "$CHROM_FILE" ]] || fail "no sequences matched -c '$CHR_REGEX' — check your header naming"
ok "$CHROM_FILE ($(wc -l < "$CHROM_FILE") chromosomes; $(( $(wc -l < "$REFERENCE.fai") - $(wc -l < "$CHROM_FILE") )) unplaced excluded)"

echo "=================================================================="
echo "4. BAM list validation"
echo "=================================================================="
[[ -s "$BAM_LIST" ]] || fail "BAM list missing/empty: $BAM_LIST"
N_BAM=$(grep -c . "$BAM_LIST")
MISSING=0
while read -r b; do [[ -z "$b" ]] && continue; [[ -f "$b" ]] || { echo "  missing: $b"; MISSING=$((MISSING+1)); }; done < "$BAM_LIST"
[[ "$MISSING" -eq 0 ]] || fail "$MISSING BAM(s) in the list do not exist"
ok "$N_BAM BAMs, all present"

xargs -a "$BAM_LIST" samtools quickcheck || fail "one or more BAMs are truncated (samtools quickcheck)"
ok "all BAMs pass quickcheck"

# @SQ vs reference: the classic '0 sites retained' trap
REF_FIRST=$(cut -f1 "$REFERENCE.fai" | head -1)
BAM_FIRST=$(samtools view -H "$(head -1 "$BAM_LIST")" | awk '$1=="@SQ"{sub("SN:","",$2); print $2; exit}')
[[ "$REF_FIRST" == "$BAM_FIRST" ]] || fail "BAM contig '$BAM_FIRST' != reference contig '$REF_FIRST' — BAMs were mapped to a different assembly"
ok "BAM @SQ names match the reference ($BAM_FIRST)"

# every BAM should come from the same assembly
BAD_ASM=$(while read -r b; do [[ -z "$b" ]] && continue
    s=$(samtools view -H "$b" | awk '$1=="@SQ"{sub("SN:","",$2); print $2; exit}')
    [[ "$s" == "$REF_FIRST" ]] || echo "$b"; done < "$BAM_LIST" | wc -l)
[[ "$BAD_ASM" -eq 0 ]] || fail "$BAD_ASM BAM(s) are from a different assembly — check for leftovers of an older mapping"
ok "all BAMs share one assembly"

echo "=================================================================="
echo "5. Clearing stale regenerable files"
echo "=================================================================="
# these are rebuilt by the driver; stale copies from an older assembly or a
# crashed run silently corrupt everything downstream.
for f in "$OUT_DIR"/metadata.tsv "$OUT_DIR"/random_10M_sites.list \
         "$OUT_DIR"/random_10M_sites.list.* "$OUT_DIR"/random_10M_contigs.txt; do
    [[ -e "$f" ]] && { rm -f "$f"; echo "  removed $(basename "$f")"; }
done
for d in pop_map_lists site_year_map_lists year_map_lists; do
    [[ -d "$OUT_DIR/$d" ]] || continue
    n=$(ls "$OUT_DIR/$d"/*.list 2>/dev/null | wc -l)
    [[ "$n" -gt 0 ]] && { rm -f "$OUT_DIR/$d"/*.list; echo "  cleared $n list(s) from $d/"; }
done
ok "stale files cleared (kept .bak/.withreps and singleton_excluded/)"

echo ""
echo "=================================================================="
echo "PREFLIGHT PASSED"
echo "  Reference : $REFERENCE"
echo "  Chromosomes: $CHROM_FILE"
echo "  BAMs      : $N_BAM"
echo ""
echo "Next:"
echo "  sbatch scripts/p06_angsd_driver.sh -o $OUT_DIR \\"
echo "    -r $REFERENCE \\"
echo "    -b $BAM_LIST \\"
echo "    -i 0.5 -P 1 -t 8 -q 20 -p 20"
echo "=================================================================="
