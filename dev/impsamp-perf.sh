#!/usr/bin/env bash
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
L=impsamp-out/perf.txt
: > $L
for cfg in "5000 5 2 2 64 1" "5000 5 1 2 64 1" "5000 5 3 3 64 0" "2000 5 2 40 64 1" "1000 5 2 2 500 0"; do
  for rep in 1 2; do
    for arm in after before; do
      REV_ARM=$arm "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" impsamp-perf.R $cfg 2>&1 | grep -E "^(RESULT|WHOLE|Error)" >> $L
    done
  done
done
echo done
