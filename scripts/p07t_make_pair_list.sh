#!/bin/bash
# Lists the BAMs of animals sequenced twice: the 10 raw/amplified 2023 pairs (o = raw DNA) and the
# 2023 animals re-extracted from tissue in the 2024 run (SI19T, SI14aT, SI29aT).
# Run on the login node from /mnt/parscratch/users/bi4og:  bash scripts/p07t_make_pair_list.sh
ALL="/mnt/parscratch/users/bi4og/physalia/angsd_results/sex_det/window_depth/bam_list.txt"   # the 195 BAMs from step A
OUT="physalia/sample_lists/bam_list_batch_pairs.list"
grep -E '/[0-9]+o?-(CAI6|NS1|NS14|NS30|NS32|NS40|NS41|SI14|SI18|SI19|SI29|SI59|SI7)(o|T|aT)?_merged' "$ALL" | sort > "$OUT"
echo "$(wc -l < "$OUT") BAMs written to $OUT (expect 26):"
xargs -a "$OUT" -n1 basename | sed 's/_merged.*//'
