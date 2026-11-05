#!/usr/bin/env bash
# After the ridge row's bound: that file on each BLAS and in the gated suite, the two
# mutants it can move, then R CMD check last. Nothing here installs.
set -u
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/ncgp-out
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
bash dev/ciharden-openblas.sh 0.3.26 lapack > $O/ob26-build.txt 2>&1
bash dev/ciharden-openblas.sh 0.3.32 lapack > $O/ob32-build.txt 2>&1
(
  export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
  F=tests/testthat/test-gp-start.R
  OPENBLAS_NUM_THREADS=4 dev/ciharden-out/Rob0.3.26-lapack/bin/Rscript.exe dev/ncgp-runtest.R frmtmb $F > $O/ob26/test-gp-start.txt 2>&1 &
  OPENBLAS_NUM_THREADS=4 dev/ciharden-out/Rob0.3.32-lapack/bin/Rscript.exe dev/ncgp-runtest.R frmtmb $F > $O/ob32/test-gp-start.txt 2>&1 &
  OPENBLAS_NUM_THREADS=1 dev/ciharden-out/Rob0.3.26-lapack/bin/Rscript.exe dev/ncgp-runtest.R frmtmb $F > $O/ob26/test-gp-start-t1.txt 2>&1 &
  OPENBLAS_NUM_THREADS=1 dev/ciharden-out/Rob0.3.32-lapack/bin/Rscript.exe dev/ncgp-runtest.R frmtmb $F > $O/ob32/test-gp-start-t1.txt 2>&1 &
  "$RS" dev/ncgp-runtest.R frmtmb $F > $O/gated/frmtmb--test-gp-start.txt 2>&1 &
  wait
)
rm -rf dev/ciharden-out
grep -h '^RESULT' $O/gated/*.txt > $O/gated.txt
"$RS" dev/ncgp-mut.R nopolish > $O/mut/nopolish.txt 2>&1 &
"$RS" dev/ncgp-mut.R control > $O/mut/control.txt 2>&1 &
wait
grep -h NCGPM $O/mut/*.txt > $O/mut.txt
echo STAGE9-FILES
