#!/bin/bash
# The Laplace simulation: 200 datasets per design, seeds 1..200.
cd /c/Users/adf44/source/r/frmtmb-wt-latres
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/review-latres-out/sim
for d in pois bern pois_ar9 pois_s1 bern_long binom3 pois_big; do
  for s in 1 51 101 151; do
    e=$((s+49))
    echo "$d $s $e dev/review-latres-out/sim/$d-$s.tsv"
  done
done | xargs -P 6 -L 1 "$R" dev/review-latres-sim.R > dev/review-latres-out/sim.log 2>&1
