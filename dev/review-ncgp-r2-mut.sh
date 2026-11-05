#!/usr/bin/env bash
# Fix round 1, item 12: the lane's mutants, run again, one per process, after the cost tables.
cd /c/Users/adf44/source/r/frmtmb-wt-ncgp
O=dev/review-ncgp-out
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
until grep -q COSTDONE $O/r2-cost.txt; do sleep 20; done
export NOT_CRAN=true
for m in control nowarm onepoint noscan noorder norowcheck nopolish; do
  "$R" dev/ncgp-mut.R $m 2>&1 | grep "^NCGPM" | cut -c1-170 >> $O/r2-mut.txt &
  while [ "$(jobs -r | wc -l)" -ge 4 ]; do sleep 5; done
done
wait
echo MUTDONE >> $O/r2-mut.txt
