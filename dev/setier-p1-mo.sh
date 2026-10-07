#!/usr/bin/env bash
# Lane setier, punch round 1: copy of dev/setier-rev-mo.sh with its own
# log directory, lane arm only: the 200-seed mo() study, ref BLAS.
W=/c/Users/adf44/source/r/frmtmb-wt-setier
cd "$W"
L=dev/setier-log/p1/mo; mkdir -p $L
REF="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
lib=C:/Users/adf44/source/r/wt-setier-lib
for c in 1 2 3 4; do
  a=$(( (c-1)*50+1 )); b=$(( c*50 ))
  "$REF" dev/setier-rev-mo.R $lib "$a:$b" $L/mo-lane-$c.tsv \
    > $L/mo-lane-$c.txt 2>&1 &
done
wait
