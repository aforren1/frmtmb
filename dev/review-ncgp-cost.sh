#!/usr/bin/env bash
# Item 8: base, lane and warm one at a time, 3 repeats each, in that order per size.
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd /c/Users/adf44/source/r/frmtmb-wt-ncgp
for n in 60 200 500 1000; do for a in base warm lane; do
  "$R" dev/review-ncgp-cost.R $a $n 3 2>&1 | grep "^RCO\|^lib\|rror"
done; done
echo COSTDONE
