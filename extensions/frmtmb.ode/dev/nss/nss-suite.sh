#!/bin/bash
# One test file per R process, which is the house rule: a whole-suite
# run in one process has repeatedly hidden state leakage here.
#
# Usage: nss-suite.sh <package-dir> [tag]
# Script path: extensions/frmtmb.ode/dev/nss/nss-suite.sh
set -u
PKG=${1:?package directory}
TAG=${2:-suite}
OUT=C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss
R="/c/Program Files/R/R-4.6.1/bin"
export R_LIBS="C:/Users/adf44/source/r/nsslib;C:/Users/adf44/source/r/rellib-0552;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true

DIR="$PKG/tests/testthat"
: > "$OUT/nss-$TAG-summary.txt"
for f in "$DIR"/test-*.R; do
  nm=$(basename "$f")
  "$R/Rscript.exe" "$OUT/nss-one.R" "$DIR" "$nm" \
      > "$OUT/nss-$TAG-$nm.log" 2>&1
  grep "^FILE" "$OUT/nss-$TAG-$nm.log" \
      >> "$OUT/nss-$TAG-summary.txt" ||
    echo "FILE $nm NO SUMMARY LINE (the process did not finish)" \
      >> "$OUT/nss-$TAG-summary.txt"
done
cat "$OUT/nss-$TAG-summary.txt"
