#!/usr/bin/env bash
# Fix round 1 check: the lane's affected test files on the lane library,
# and the changed ones on rellib-r13, one file per process, 5 at a time.
cd /c/Users/adf44/source/r/frmtmb-wt-mocode || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
out=dev/review-mocode-log/fix1-files
mkdir -p "$out"
{
  for f in tests/testthat/test-se-check.R tests/testthat/test-mo-coding.R \
           tests/testthat/test-emmeans.R tests/testthat/test-ce-*.R; do
    echo "lane frmtmb $f"
  done
  for f in extensions/frmtmb.sample/tests/testthat/test-mo-coding.R \
           extensions/frmtmb.sample/tests/testthat/test-conditional-effects-draws.R; do
    echo "lane frmtmb.sample $f"
  done
  for f in tests/testthat/test-se-check.R tests/testthat/test-mo-coding.R \
           tests/testthat/test-emmeans.R; do
    echo "base frmtmb $f"
  done
  echo "base frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-mo-coding.R"
} | xargs -P 5 -L 1 bash -c '"$0" dev/mocode-run-test.R "$1" "$2" "$3" > "'"$out"'/$1-$2-$(basename "$3" .R).txt" 2>&1' "$R"
grep -H '^RESULT' "$out"/*.txt | sed 's#.*/##'
