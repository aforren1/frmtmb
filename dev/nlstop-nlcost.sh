#!/usr/bin/env bash
# Lane nlstop: wall time of the nonlinear test files on rellib-r20 and on
# the lane's build, alternating base and lane per file, 3 rounds, one R
# process at a time. Prints one TIME line per run; the summary takes the
# minimum per file and build.
# Usage: bash dev/nlstop-nlcost.sh > dev/nlstop-out/nlcost.txt
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
files="test-nl-start.R test-nl.R test-nl-family.R test-nlf.R"
for round in 1 2 3; do
  for f in $files; do
    [ -f "tests/testthat/$f" ] || continue
    for arm in base lane; do
      if [ "$arm" = lane ]; then
        export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlstop-lib
      else
        unset FRMTMB_RELLIB
      fi
      t0=$(date +%s.%N)
      r=$("$RS" dev/release/run-tests.R frmtmb "tests/testthat/$f" 2>&1 |
            grep -E "^RESULT|^lib")
      t1=$(date +%s.%N)
      echo "TIME round $round $f $arm $(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b - a}') | $(echo "$r" | tr '\n' ' ')"
    done
  done
done
