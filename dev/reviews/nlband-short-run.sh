#!/usr/bin/env bash
# Short check of lane nlband, fix round 1: one arm per call.
# bash dev/reviews/nlband-short-run.sh <arm> <sections> <log> [check|extra]
cd /c/Users/adf44/source/r/frmtmb-wt-nlband || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true R_PARALLEL_PORT=$((11500 + RANDOM % 400))
NLBAND_ARM=$1 "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" \
  dev/reviews/nlband-short-${4:-check}.R "$2" > dev/reviews/nlband-short-out/$3 2>&1
echo "EXIT $?" >> dev/reviews/nlband-short-out/$3
