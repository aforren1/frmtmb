#!/bin/bash
# Lane nlstop2: run one R script of the lane from the worktree, log to
# dev/nlstop2-out/<log>. Usage: NLSTOP2_ARM=before|after
#   bash dev/nlstop2-run.sh <script.R> <log> [args...]
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop2 || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
s="$1"; log="$2"; shift 2
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$s" "$@" \
  > "dev/nlstop2-out/$log" 2>&1
echo "EXIT $? $s $log" >> "dev/nlstop2-out/$log"
