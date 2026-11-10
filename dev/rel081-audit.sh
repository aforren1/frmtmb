#!/usr/bin/env bash
# The translator's per-row audit beside the 0.81.0 record (the
# translator review's procedure note: `eq_degenerate` holds the build's
# estimates, so a build that moves an estimate onto a boundary needs the
# audit again). Pairs: every pair whose recorded level holds L1, plus
# p069 (L1 on its Laplace fit), as the lane chose them; compared with
# the lane's list. 2 jobs (one R process per pair).
#
#   bash dev/rel081-audit.sh
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/rel081-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel081 || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export FRMTMB_EQUIV_LIB="$L/rellib-r20;$L/rellib-r19"
O=dev/rel081-audit-out
G=dev/rel081-log
ids=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) { h = $i; gsub(/"/, "", h); c[h] = i } ; next }
  { id = $c["id"]; lv = $c["level"]; gsub(/"/, "", id); gsub(/"/, "", lv)
    if (lv ~ /L1/ || id == "p069") print id }' dev/equiv-data/results.tsv | sort -u)
printf '%s\n' $ids > $G/audit-ids.txt
echo "pairs: $(grep -c . $G/audit-ids.txt)"
tr ' ' '\n' < $L/frmtmb-wt-translator/dev/translator-audit-ids.txt | grep . | sort -u > $G/audit-ids-lane.txt
echo "lane's list: $(grep -c . $G/audit-ids-lane.txt) (+ p069); only now: $(comm -23 $G/audit-ids.txt $G/audit-ids-lane.txt | tr '\n' ' '); only the lane's: $(comm -13 $G/audit-ids.txt $G/audit-ids-lane.txt | tr '\n' ' ')"
rm -rf $O
bash dev/translator-audit.sh $O 2 $ids
"$RS" dev/translator-audit-summary.R $O > $G/audit-summary.txt 2>&1
cat $G/audit-summary.txt
grep -h '^lib' $O/logs/*.log | sort | uniq -c
echo "AUDIT END"
