#!/bin/sh
# Run one review script from the worktree root: run.sh <log> <script> [args]
cd C:/Users/adf44/source/r/frmtmb-wt-nlstop
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true R_PARALLEL_PORT=$((11000 + $$ % 1000))
log=$1; shift
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$@" > dev/reviews/nlstop-short-out/$log 2>&1
echo "EXIT $?" >> dev/reviews/nlstop-short-out/$log
