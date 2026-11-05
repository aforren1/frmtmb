#!/bin/bash
cd /c/Users/adf44/source/r/frmtmb-wt-latres
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/review-latres-out/m1b
for d in zip zip_hi zip_lo hp zib cum3 cum5 cum7 nb_m1 nb_m5 nb_m20 nb_s1 beta_lo binom4 ln_s2 ln_s05 sn_u1 sn_u10; do
  echo "$d 3001 3040 dev/review-latres-out/m1b/$d.tsv"
done | xargs -P 3 -L 1 "$R" dev/review-latres-m1b.R > dev/review-latres-out/m1b.log 2>&1
