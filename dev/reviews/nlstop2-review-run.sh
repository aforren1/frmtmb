#!/bin/bash
# Review of lane nlstop2: run one review script, log to
# dev/reviews/nlstop2-review-out/<log>. Usage:
#   REV_ARM=before|after bash dev/reviews/nlstop2-review-run.sh <s.R> <log> [args]
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop2 || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
mkdir -p dev/reviews/nlstop2-review-out
s="$1"; log="$2"; shift 2
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$s" "$@" \
  > "dev/reviews/nlstop2-review-out/$log" 2>&1
echo "EXIT $? $s $log" >> "dev/reviews/nlstop2-review-out/$log"
