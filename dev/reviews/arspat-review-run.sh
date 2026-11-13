#!/bin/bash
# usage: bash dev/reviews/arspat-review-run.sh <script.R> <log> [args]
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-arspat/dev/stan-cache
cd C:/Users/adf44/source/r/frmtmb-wt-arspat
s=$1; l=$2; shift 2
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$s" "$@" > "$l" 2>&1
echo "exit $?" >> "$l"
