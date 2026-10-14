#!/bin/sh
# Fix round 1 provenance run, against the final installs of core, eam and
# sample: the full eight-package suite, the gated tier, the drmTMB gate,
# the family test files under OpenBLAS 0.3.26 with LAPACK. Each step's
# summary goes to dev/famuni-out/final-*.txt.
WT=/c/Users/adf44/source/r/frmtmb-wt-famuni
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd "$WT" || exit 1
ls -l --time-style=full-iso C:/Users/adf44/source/r/wt-famuni-lib/frmtmb/DESCRIPTION \
  C:/Users/adf44/source/r/wt-famuni-lib/frmtmb.eam/DESCRIPTION \
  C:/Users/adf44/source/r/wt-famuni-lib/frmtmb.sample/DESCRIPTION \
  > dev/famuni-out/final-installs.txt
sh dev/famuni-suite.sh full 6 > dev/famuni-out/final-suite.driverlog 2>&1
"$RS" dev/famuni-summary.R dev/famuni-suite-full > dev/famuni-out/final-suite.txt 2>&1
sh dev/famuni-gated.sh gated 6 > dev/famuni-out/final-gated.driverlog 2>&1
"$RS" dev/famuni-summary.R dev/famuni-suite-gated > dev/famuni-out/final-gated.txt 2>&1
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$WT/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
( cd tests/testthat && "$RS" ../../dev/famuni-test1.R frmtmb test-drmtmb-agreement.R \
    > ../../dev/famuni-out/final-drmtmb.txt 2>&1 )
bash dev/ciharden-openblas.sh 0.3.26 lapack > dev/famuni-out/final-openblas-build.txt 2>&1
OB="$WT/dev/ciharden-out/Rob0.3.26-lapack/bin/Rscript.exe"
export OPENBLAS_NUM_THREADS=4
: > dev/famuni-out/final-openblas.txt
for f in test-families-extreme.R test-simulate-density.R; do
  ( cd tests/testthat && "$OB" ../../dev/famuni-test1.R frmtmb $f > ../../dev/famuni-out/final-ob-$f.txt 2>&1 )
  grep RESULT dev/famuni-out/final-ob-$f.txt >> dev/famuni-out/final-openblas.txt
done
for f in test-gev-tmp-xi.R test-default-priors-brms.R; do
  ( cd extensions/frmtmb.sample/tests/testthat && "$OB" ../../../../dev/famuni-test1.R frmtmb.sample $f > ../../../../dev/famuni-out/final-ob-$f.txt 2>&1 )
  grep RESULT dev/famuni-out/final-ob-$f.txt >> dev/famuni-out/final-openblas.txt
done
( cd extensions/frmtmb.eam/tests/testthat && "$OB" ../../../../dev/famuni-test1.R frmtmb.eam test-wiener-links.R > ../../../../dev/famuni-out/final-ob-test-wiener-links.R.txt 2>&1 )
grep RESULT dev/famuni-out/final-ob-test-wiener-links.R.txt >> dev/famuni-out/final-openblas.txt
rm -rf dev/ciharden-out
echo FINAL DONE > dev/famuni-out/final-done.txt
