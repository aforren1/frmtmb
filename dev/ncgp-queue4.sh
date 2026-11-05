#!/usr/bin/env bash
# After the draws-order fix in frmtmb.sample: that package's gated tier and its
# three one-session runs, its two OpenBLAS files, and R CMD check of both packages.
# The sampler's own read on the design at 1e-2 runs beside this.
set -u
rm -f dev/ncgp-out/gated/frmtmb.sample--*.txt
FRMTMB_BRMS_FIT_TESTS=true bash dev/ncgp-par.sh dev/ncgp-out/gated frmtmb.sample
echo GATED-DONE
(
  export NOT_CRAN=true
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-ncgp-lib
  export TMP="${TMP:-$TEMP}"
  bash dev/release/run-sessions.sh ncgpsess2 3 frmtmb.sample > dev/ncgp-out/sess2-driver.txt 2>&1
  cp dev/release/ncgpsess2.log dev/ncgp-out/ncgpsess2.log 2>/dev/null
)
echo SESS-DONE
bash dev/ciharden-openblas.sh 0.3.26 lapack > dev/ncgp-out/ob26-build.txt 2>&1
bash dev/ciharden-openblas.sh 0.3.32 lapack > dev/ncgp-out/ob32-build.txt 2>&1
(
  rt=/c/rtools45
  export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
  export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true OPENBLAS_NUM_THREADS=4
  export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
  export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
  for v in 26 32; do
    RB=dev/ciharden-out/Rob0.3.$v-lapack/bin/Rscript.exe
    for f in test-gp-by-draws test-reparam; do
      "$RB" dev/ncgp-runtest.R frmtmb.sample \
        "extensions/frmtmb.sample/tests/testthat/$f.R" \
        > "dev/ncgp-out/ob$v/$f.txt" 2>&1 &
    done
    wait
  done
)
echo OB-DONE
bash dev/ncgp-check.sh extensions/frmtmb.sample > dev/ncgp-out/check-sample.txt 2>&1
rm -rf dev/ciharden-out
bash dev/ncgp-check.sh . > dev/ncgp-out/check-core.txt 2>&1
echo CHECK-DONE
