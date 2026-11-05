#!/usr/bin/env bash
# Short check of fix round 1: my review's scripts again, on the fixed build.
O=dev/review-refuse-out/f1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
{
for f in student exgaussian skew_normal gen_extreme_value beta Gamma inverse.gaussian von_mises; do echo "decision $f"; done
for f in gaussian student beta Gamma exponential exgaussian von_mises; do echo "mi-exact $f"; done
for s in brms-exg misc misc2 misc3 misc4 mm-weights mm2 mm-hand const mi-edges mi-methods; do echo "$s x"; done
for c in Gamma von_mises; do echo "sample $c"; done
} | xargs -P 5 -L1 bash -c '"$0" dev/review-refuse-$1.R $2 > dev/review-refuse-out/f1/$1-$2.log 2>&1' "$R"
