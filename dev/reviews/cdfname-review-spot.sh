#!/bin/bash
# Review spot run: eam and survival family tests, two R processes at once.
cd /c/Users/adf44/source/r/frmtmb-wt-cdfname
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
run() { "$R" dev/cdfname-run-test.R "$1" "extensions/$1/tests/testthat/$2" > "dev/reviews/cdfname-review-spot-$2.log" 2>&1; }
run frmtmb.eam test-rdm-gng.R & run frmtmb.survival test-rp-floored.R; wait
run frmtmb.eam test-wiener-cdf.R & run frmtmb.survival test-royston-parmar.R; wait
run frmtmb.eam test-family.R & run frmtmb.survival test-saved-fits.R; wait
grep -h 'lib:\|RESULT\|FAILED' dev/reviews/cdfname-review-spot-*.log
