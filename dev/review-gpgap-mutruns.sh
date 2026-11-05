#!/usr/bin/env bash
# item 6: 4 of the lane's 17 mutants, from a copy of its script on my own library.
set -u
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true
for m in control thr_double q_med no_trycatch first_only; do
  "$RS" dev/review-gpgap-mutrun.R $m > dev/review-gpgap-out/f1-mut-$m.txt 2>&1 &
done
wait
