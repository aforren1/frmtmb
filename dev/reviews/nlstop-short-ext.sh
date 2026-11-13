#!/bin/sh
# Short check: the extension test files that set a prior on class Intercept,
# lane core ahead of rellib-r20, one file per process, 3 at a time.
cd C:/Users/adf44/source/r/frmtmb-wt-nlstop
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlstop-lib
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
out=dev/reviews/nlstop-short-out/ext.txt; : > $out
grep -rl 'class = "Intercept"' extensions/*/tests/testthat/*.R | while read f; do
  p=$(echo $f | cut -d/ -f2)
  echo "$p $f"
done | xargs -P 3 -L 1 sh -c 'R_PARALLEL_PORT=$((13000 + $$ % 1000)) "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/release/run-tests.R $0 $1 2>&1 | grep -E "^lib:|^RESULT|^DETAIL" >> dev/reviews/nlstop-short-out/ext.txt'
echo DONE >> $out
