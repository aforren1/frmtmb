#!/usr/bin/env bash
# The OpenBLAS runs (r078-ob26, r078-ob32: dev/rel078-log/ob.jobs under
# OpenBLAS 0.3.26 and 0.3.32 with LAPACK, 4 threads) against the
# reference BLAS run of the same files in the gated tier (r078-ref): the
# RESULT line of each file, and the blas: line of each OpenBLAS log.
#   bash dev/rel078-ob-compare.sh > dev/rel078-log/ob-compare.txt
cd /c/Users/adf44/source/r/frmtmb-wt-rel078 || exit 1
L=dev/ciharden-log
n=0; d=0
while read -r p f; do
  b=$(basename "$f" .R)
  ref=$(grep -h "^RESULT" "$L/r078-ref-suite/$p--$b.txt" | tail -1)
  for v in ob26 ob32; do
    o=$(grep -h "^RESULT" "$L/r078-$v-suite/$p--$b.txt" | tail -1)
    n=$((n + 1))
    if [ "$o" != "$ref" ]; then d=$((d + 1)); echo "DIFFERS $v $p $b: $o | ref: $ref"; fi
  done
done < dev/rel078-log/ob.jobs
for v in ob26 ob32; do
  echo "$v: blas lines naming the OpenBLAS copy: $(grep -h '^blas:' $L/r078-$v-suite/*.txt | grep -c 'ciharden-out/Rob0.3.')" \
       "of $(ls $L/r078-$v-suite/*.txt | wc -l)"
done
echo "compared $n file runs with the reference run, differ: $d"
