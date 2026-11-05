#!/usr/bin/env bash
# After dev/nlfix-final1.sh: the one-session runs, the OpenBLAS pass on the new
# test files, and R CMD check. Logs under dev/nlfix-*.
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
LIB=C:/Users/adf44/source/r/wt-nlfix-lib
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export NOT_CRAN=true
# rerun what changed since the chain read it
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" -e 'invisible(0)'
export FRMTMB_RELLIB="$LIB"
bash dev/release/run-sessions.sh nlfixsess 8 > dev/nlfix-sessions.out 2>&1
cp dev/release/nlfixsess.log dev/nlfix-sessions.log
rm -rf dev/release/nlfixsess-sessions.keep
unset FRMTMB_RELLIB
# OpenBLAS with LAPACK, as the runners have it, on the new test files
bash dev/ciharden-openblas.sh 0.3.26 lapack > dev/nlfix-openblas-build.log 2>&1
RB="dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
for f in tests/testthat/test-nl-family.R tests/testthat/test-nl-stan-funs.R; do
  OPENBLAS_NUM_THREADS=4 "$RB" dev/nlfix-testrun.R frmtmb "$f" dev/nlfix-log/ob-$(basename "$f" .R).txt > /dev/null 2>&1
done
FRMTMB_BRMS_FIT_TESTS=true OPENBLAS_NUM_THREADS=4 "$RB" dev/nlfix-testrun.R frmtmb tests/testthat/test-nl-brms.R dev/nlfix-log/ob-test-nl-brms.txt > /dev/null 2>&1
rm -rf dev/ciharden-out
# R CMD check once per changed package
bash dev/nlfix-check.sh . > dev/nlfix-check-core.out 2>&1 &
bash dev/nlfix-check.sh extensions/frmtmb.sample > dev/nlfix-check-sample.out 2>&1 &
wait
echo FINAL2 DONE > dev/nlfix-final2.done
