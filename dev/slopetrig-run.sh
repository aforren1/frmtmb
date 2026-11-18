#!/usr/bin/env bash
# Lane slopetrig: run one R script with the lane's environment.
# Usage: bash dev/slopetrig-run.sh <log> <script> [args...]
log="$1"; shift
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-slopetrig/dev/stan-cache
cd /c/Users/adf44/source/r/frmtmb-wt-slopetrig || exit 1
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$@" > "$log" 2>&1
echo "EXIT $?" >> "$log"
