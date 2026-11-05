#!/bin/bash
cd /c/Users/adf44/source/r/frmtmb-wt-latres
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/review-latres-out/sim2
{
for d in gamma_T nb_T; do for s in 1 26 51 76; do e=$((s+24)); echo "$d $s $e dev/review-latres-out/sim2/$d-$s.tsv"; done; done
for d in gamma_F probit_F pois_F bern_F; do for s in 1 11 21; do e=$((s+9)); echo "$d $s $e dev/review-latres-out/sim2/$d-$s.tsv"; done; done
} | xargs -P 5 -L 1 "$R" dev/review-latres-sim2.R > dev/review-latres-out/sim2.log 2>&1
