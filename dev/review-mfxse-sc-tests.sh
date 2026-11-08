#!/bin/sh
# Review of lane mfxse, short check: the five test files on the fix build,
# process, one process at a time, on the lane library.
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
W=/c/Users/adf44/source/r/frmtmb-wt-mfxse
out=$W/dev/review-mfxse-out/sc/tests
mkdir -p "$out"
for f in test-marginaleffects-re.R test-marginaleffects-args.R \
         test-me-newdata.R test-v07.R test-compat.R; do
  "$R" $W/dev/review-mfxse-run-test.R $W/tests/testthat/$f lane \
    > "$out/${f%.R}.log" 2>&1
done
grep -h -E "^(lib:|RESULT)" "$out"/*.log
