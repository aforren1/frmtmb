#!/usr/bin/env bash
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
mkdir -p review-impsamp-out/unbiased-log
for d in smallpois nearsep bigsd_bern t_outlier wtbin re_sigma probit nb; do echo "$d"; done |
  xargs -P 4 -L 1 bash -c '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" review-impsamp-unbiased.R "$0" 200 200 > review-impsamp-out/unbiased-log/"$0".log 2>&1'
echo done
