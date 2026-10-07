#!/usr/bin/env bash
# Lane setier: alternate the base and lane arms of dev/setier-cost.R,
# five rounds, sequentially, into dev/setier-log/cost.txt.
cd /c/Users/adf44/source/r/frmtmb-wt-setier
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r/wt-setier-lib
: > dev/setier-log/cost.txt
for r in 1 2 3 4 5; do
  "$RS" dev/setier-cost.R base $r >> dev/setier-log/cost.txt 2>&1
  "$RS" dev/setier-cost.R $L $r >> dev/setier-log/cost.txt 2>&1
done
