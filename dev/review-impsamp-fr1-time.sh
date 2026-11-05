#!/usr/bin/env bash
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
L=review-impsamp-out/fr1/time.txt
: > $L
for rep in 1 2 3 4; do
  for arm in after before; do
    REV_ARM=$arm "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" review-impsamp-fr1-time.R 2>&1 | grep "^TIME" >> $L
  done
done
echo done
