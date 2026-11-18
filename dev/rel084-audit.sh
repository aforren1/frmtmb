#!/usr/bin/env bash
# The translator audit beside the recorded 0.84.0 survey, as at 0.83.0
# (dev/rel083-audit.sh): every pair the record reads at L1, and p069, on
# this build (eq_degenerate is build-specific), 2 jobs.
#
#   bash dev/rel084-audit.sh
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel084 || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export FRMTMB_EQUIV_LIB="$L/rellib-r23;$L/rellib-r22"
O=dev/rel084-audit-out
G=dev/rel084-log
ids=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) { h = $i; gsub(/"/, "", h); c[h] = i } ; next }
  { id = $c["id"]; lv = $c["level"]; gsub(/"/, "", id); gsub(/"/, "", lv)
    if (lv ~ /L1/ || id == "p069") print id }' dev/equiv-data/results.tsv | sort -u)
printf '%s\n' $ids > $G/audit-ids.txt
echo "pairs: $(grep -c . $G/audit-ids.txt)"
rm -rf $O
bash dev/translator-audit.sh $O 2 $ids
"$RS" dev/translator-audit-summary.R $O > $G/audit-summary.txt 2>&1
cat $G/audit-summary.txt
grep -h '^lib' $O/logs/*.log | sort | uniq -c
echo "AUDIT END"
