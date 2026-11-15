#!/usr/bin/env bash
# Lane gpfield: run one R script with this lane's environment, output to
# dev/gpfield-out/<log>. Usage: gpfield-run.sh <log> <script> [args...]
WT=C:/Users/adf44/source/r/frmtmb-wt-gpfield
export PATH="/c/Program Files/R/R-4.6.1/bin:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
log=$1; shift
mkdir -p $WT/dev/gpfield-out
Rscript "$@" > $WT/dev/gpfield-out/$log 2>&1
echo "EXIT $? $log"
