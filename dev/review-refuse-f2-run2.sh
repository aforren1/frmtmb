#!/usr/bin/env bash
{
for f in Gamma exponential beta exgaussian student lognormal von_mises gaussian; do echo "table1 $f"; done
for s in 0 0.02 0.05 0.1 0.2 0.4 0.7; do echo "table2 $s"; done
echo "vm x"
} | xargs -P 7 -L1 bash -c '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/review-refuse-f2-rule.R $0 $1 > dev/review-refuse-out/f2/f2-rule-$0-$1.log 2>&1'
touch dev/review-refuse-out/f2/DONE2
