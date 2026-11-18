#!/usr/bin/env bash
# Review of lane slopetrig: run one R script from the worktree with the
# lane environment. Usage: slopetrig-review-run.sh <port> <log> <script> [args]
cd /c/Users/adf44/source/r/frmtmb-wt-slopetrig || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_PARALLEL_PORT="$1"
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-slopetrig/dev/stan-cache"
log="$2"
shift 2
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$@" > "$log" 2>&1
echo "EXIT $?" >> "$log"
