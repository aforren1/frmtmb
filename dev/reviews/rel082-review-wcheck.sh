#!/usr/bin/env bash
# Release review of 0.82.0: lane wcheck's 58-case table on rellib-r21,
# one R process at a time, diffed against the lane's fr3 table.
set -u
cd C:/Users/adf44/source/r/frmtmb-wt-rel082
export PATH="/c/Program Files/R/R-4.6.1/bin:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
O=dev/reviews/rel082-review-out/wcheck
W=C:/Users/adf44/source/r/frmtmb-wt-wcheck/dev
mkdir -p $O
for p in mv cens prior; do
  Rscript.exe --vanilla dev/wcheck-cases.R C:/Users/adf44/source/r/rellib-r21 $p > $O/cases-$p-lane.txt 2>&1
  cp $W/wcheck-out/cases-$p-base.txt $O/
done
sed "s#dev/wcheck-out/#$O/#g" dev/wcheck-table.R > $O/wcheck-table.R
Rscript.exe --vanilla $O/wcheck-table.R > $O/table.txt 2>&1
diff $W/wcheck-out/fr3/table.txt $O/table.txt > $O/table.diff
echo "diff exit $?" >> $O/table.diff
grep -h "^lib:" $O/cases-*-lane.txt > $O/libs.txt
echo END > $O/done
