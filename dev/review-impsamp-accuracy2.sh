#!/usr/bin/env bash
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
for d in zip hurdle betafam cens_gauss trunc_pois scale1e5 pois1e5 year tiny_x_slope huge_group zero_wt; do for s in 1 2 3; do echo "$d $s"; done; done |
  xargs -P 3 -L 1 bash -c '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" review-impsamp-accuracy2.R "$0" "$1" > review-impsamp-out/accuracy2-log/"$0"-"$1".log 2>&1'
echo done
