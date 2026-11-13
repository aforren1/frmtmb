#!/bin/sh
# Usage: silent-review-run.sh <logname> <script> args...
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
log=$1; shift
cd /c/Users/adf44/source/r/frmtmb-wt-silent
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$@" > dev/reviews/silent-review-log/$log.log 2>&1
echo "EXIT $?" >> dev/reviews/silent-review-log/$log.log
