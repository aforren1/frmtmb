#!/bin/sh
# One test file per R process, sequentially. Reviewer's batch.
# Usage: sh dev/asrev-batch.sh <logfile> <pkgdir> <file> [<file> ...]
R="C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
LOG="$1"; shift
PKG="$1"; shift
cd "C:/Users/adf44/source/r/frmtmb-wt-argspell" || exit 1
: > "$LOG"
for f in "$@"; do
  echo "######## $PKG $f" >> "$LOG"
  "$R" dev/asrev-run1.R "$PKG" "^$f\$" >> "$LOG" 2>&1
done
grep -E "^(########|ASREV)" "$LOG"
