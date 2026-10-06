#!/bin/sh
# Reviewer test driver: one R process per file, lane arm gated, plus
# the seen-to-fail arm on rellib-r5 for three files.
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
W=/c/Users/adf44/source/r/frmtmb-wt-gpby
cd $W
mkdir -p dev/gpby-rev-tests
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
run() { # pkg dir file arm
  out="dev/gpby-rev-tests/$4-$1-$(basename $3 .R).txt"
  "$R" dev/gpby-rev-runtest.R $1 $2/$3 $4 > $out 2>&1
}
C=tests/testthat
S=extensions/frmtmb.sample/tests/testthat
P=extensions/frmtmb.spline/tests/testthat
for f in test-gp-by.R test-fd-chain.R test-gp-multidim.R test-emmeans.R \
         test-predict-re-uncertainty.R test-prior-compat.R \
         test-brms-agreement.R test-get-prior-route.R test-compat.R \
         test-brms-priors.R test-ce-levels.R test-ce-parity.R test-brms-likelihood.R \n         test-brms-suite-methods.R test-brms-suite-emmeans.R test-perf.R; do
  [ -f $C/$f ] && run frmtmb $C $f lane &
done
wait
for f in test-gp-by-draws.R test-prior-route.R test-gp.R; do
  [ -f $S/$f ] && run frmtmb.sample $S $f lane &
done
for f in test-curve.R test-difference.R test-new-levels.R test-deriv.R; do
  [ -f $P/$f ] && run frmtmb.spline $P $f lane &
done
run frmtmb.eam extensions/frmtmb.eam/tests/testthat test-family.R lane &
run frmtmb.ode extensions/frmtmb.ode/tests/testthat test-compat.R lane &
run frmtmb.latent extensions/frmtmb.latent/tests/testthat test-prior-dpar.R lane &
run frmtmb.learn extensions/frmtmb.learn/tests/testthat test-families.R lane &
run frmtmb.coupling extensions/frmtmb.coupling/tests/testthat test-conditions.R lane &
wait
run frmtmb $C test-gp-by.R base &
run frmtmb $C test-fd-chain.R base &
run frmtmb.spline $P test-difference.R base &
wait
echo ALLDONE
