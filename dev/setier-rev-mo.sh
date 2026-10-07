#!/usr/bin/env bash
# Reviewer of lane setier: the 200-seed mo() study, lane and base, ref BLAS.
W=/c/Users/adf44/source/r/frmtmb-wt-setier
cd "$W"
L=dev/setier-rev-log/mo; mkdir -p $L
REF="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for arm in lane base; do
  if [ $arm = lane ]; then lib=C:/Users/adf44/source/r/wt-setier-lib; else lib=C:/Users/adf44/source/r/rellib-r6; fi
  for c in 1 2 3 4; do
    a=$(( (c-1)*50+1 )); b=$(( c*50 ))
    "$REF" dev/setier-rev-mo.R $lib "$a:$b" $L/mo-$arm-$c.tsv > $L/mo-$arm-$c.txt 2>&1 &
  done
done
wait
