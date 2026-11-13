#!/usr/bin/env bash
# Runs one R script of the 0.82.0 release review on rellib-r21; log beside it.
# usage: rel082-review-r.sh <script.R> <log> [args...]
set -u
export PATH="/c/Program Files/R/R-4.6.1/bin:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_LIBS_USER=
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-rel082/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
s=$1; log=$2; shift 2
cd C:/Users/adf44/source/r/frmtmb-wt-rel082
Rscript.exe --vanilla "$s" "$@" > "$log" 2>&1
echo "exit $?" >> "$log"
