#!/usr/bin/env bash
# Reviewer: every RESULT line of the reviewer's test runs, per tag, with
# the package, generated from the logs. Output:
# dev/aterms2-rev-log-11-counts.txt
cd /c/Users/adf44/source/r/frmtmb-wt-aterms2/dev
for d in aterms2-rev-suite-*; do
  tag=${d#aterms2-rev-suite-}
  for f in $d/*.txt; do
    pkg=$(grep -m1 "^ARM" "$f" | sed 's#.*/\(frmtmb[.a-z]*\) *$#\1#')
    res=$(grep -m1 "^RESULT" "$f" || echo "RESULT MISSING $(basename $f)")
    echo "$tag $pkg $res"
  done
done | sort -k1,1 -k2,2 -k5,5 -k4,4
echo "files: $(ls aterms2-rev-suite-*/*.txt | wc -l)  with RESULT: $(grep -l '^RESULT' aterms2-rev-suite-*/*.txt | wc -l)"
