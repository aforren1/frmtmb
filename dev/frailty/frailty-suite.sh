#!/bin/sh
# The whole frmtmb.spline suite, ONE TEST FILE PER R PROCESS.
cd "$(dirname "$0")" || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
DIR=/c/Users/adf44/source/r/frmtmb-wt-frailty/extensions/frmtmb.spline/tests/testthat
for f in $(ls "$DIR" | grep '^test-.*\.R$'); do
  out=$("$R" frailty-runtest.R "$f" 2>&1)
  line=$(echo "$out" | grep '^FILE ')
  if [ -z "$line" ]; then
    echo "FILE $f ABORTED"
    echo "$out" | tail -6
  else
    echo "$line"
  fi
done
echo "SUITE DONE"
