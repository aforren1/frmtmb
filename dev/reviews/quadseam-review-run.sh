#!/usr/bin/env bash
# bash quadseam-review-run.sh <script> <logname> args...
set -u
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
cd /c/Users/adf44/source/r/frmtmb-wt-quadseam/dev/reviews
s=$1; l=$2; shift 2
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$s" "$@" > "quadseam-review-log/$l.txt" 2>&1
echo "exit $? $l"
