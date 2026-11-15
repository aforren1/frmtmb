#!/usr/bin/env bash
# Short check of lane gpfield: run one R script, output to dev/reviews/gpfield-short-out/<log>
WT=C:/Users/adf44/source/r/frmtmb-wt-gpfield
export PATH="/c/Program Files/R/R-4.6.1/bin:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
log=$1; shift
Rscript "$@" > $WT/dev/reviews/gpfield-short-out/$log 2>&1
echo "EXIT $? $log"
