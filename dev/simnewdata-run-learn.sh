#!/bin/sh
# The frmtmb.learn part of the ungated suite, re-run after this lane
# edited frmtmb.learn (its refusal message and one test), with the
# release runner's copy and the same environment as run-suite.ps1.
#   sh dev/simnewdata-run-learn.sh
ROOT=C:/Users/adf44/source/r/frmtmb-wt-simnewdata
cd "$ROOT" || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LOG=dev/simnewdata-log/suite-learn.log
export NOT_CRAN=true FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
: > $LOG
echo "== frmtmb.learn files ==" >> $LOG
n=0
for f in extensions/frmtmb.learn/tests/testthat/test-*.R; do
  out=$("$R" dev/simnewdata-run-tests.R frmtmb.learn "$ROOT/$f" 2>&1 | grep "^RESULT")
  if [ -n "$out" ]; then n=$((n + 1)); echo "$out" >> $LOG; else echo "NO RESULT LINE for $f" >> $LOG; fi
done
echo "LEARN ran $n" >> $LOG
