#!/bin/sh
# The family test files under OpenBLAS 0.3.26 with LAPACK routed to it
# (dev/ciharden-openblas.sh 0.3.26 lapack must have built the R copy).
# Output: dev/famuni-out/final-ob-*.txt and the RESULT lines in
# dev/famuni-out/final-openblas.txt.
WT=/c/Users/adf44/source/r/frmtmb-wt-famuni
cd "$WT" || exit 1
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$WT/dev/stan-cache"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
OB="$WT/dev/ciharden-out/Rob0.3.26-lapack/bin/Rscript.exe"
[ -x "$OB" ] || { echo "no OpenBLAS R at $OB"; exit 1; }
# (the BLAS and LAPACK routing is in final-openblas-build.txt)
run() {
  ( cd "$1" && "$OB" "$WT/dev/famuni-test1.R" "$2" "$3" \
      > "$WT/dev/famuni-out/final-ob-$3.txt" 2>&1 )
  grep RESULT "dev/famuni-out/final-ob-$3.txt" >> dev/famuni-out/final-openblas.txt
}
run tests/testthat frmtmb test-families-extreme.R
run tests/testthat frmtmb test-simulate-density.R
run extensions/frmtmb.sample/tests/testthat frmtmb.sample test-gev-tmp-xi.R
run extensions/frmtmb.sample/tests/testthat frmtmb.sample test-default-priors-brms.R
run extensions/frmtmb.eam/tests/testthat frmtmb.eam test-wiener-links.R
