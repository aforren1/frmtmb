#!/usr/bin/env bash
# Driver: every design at seeds 1..5, at most $1 R processes at once.
P=${1:-6}
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
LOGD=review-impsamp-out/accuracy-log
mkdir -p "$LOGD"
DESIGNS="unbal_bern bigsd_bern nearsep smallpois offpois wtbin probit cloglog gaussloglink distreg re_sigma re_sigma_small t_nu5 t_outlier gauss_sigma_fixed gauss_se nb"
for d in $DESIGNS; do for s in 1 2 3 4 5; do echo "$d $s"; done; done |
  xargs -P "$P" -L 1 bash -c '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" review-impsamp-accuracy.R "$0" "$1" > review-impsamp-out/accuracy-log/"$0"-"$1".log 2>&1'
echo done
