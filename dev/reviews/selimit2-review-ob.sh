#!/usr/bin/env bash
# Review of lane selimit2, item 9: the lane's SE and warning test files
# under OpenBLAS 0.3.32 (BLAS and LAPACK) with the Haswell and SandyBridge
# kernels, through the lane's dev/release/run-tests.R, two at a time.
#   bash selimit2-review-ob.sh <R copy bin dir> <out dir>
set -u
RB="$1"; OUT="$2"
L=/c/Users/adf44/source/r/frmtmb-wt-selimit2
mkdir -p "$OUT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-selimit2-lib
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$L/dev/stan-cache" OPENBLAS_NUM_THREADS=4 OPENBLAS_VERBOSE=2
jobs_list="frmtmb:tests/testthat/test-se-limit.R frmtmb:tests/testthat/test-se-tier.R frmtmb:tests/testthat/test-mi-families.R frmtmb:tests/testthat/test-se-reparam-reml.R frmtmb.eam:extensions/frmtmb.eam/tests/testthat/test-defects.R"
for ct in Haswell SandyBridge; do
  for j in $jobs_list; do
    p=${j%%:*}; f=${j#*:}
    while [ "$(jobs -rp | wc -l)" -ge 2 ]; do sleep 2; done
    OPENBLAS_CORETYPE=$ct "$RB/Rscript.exe" -e 'cat("blas:", extSoftVersion()[["BLAS"]], La_library(), "\n")' > "$OUT/$ct--$p--$(basename "$f" .R).txt" 2>&1
    OPENBLAS_CORETYPE=$ct "$RB/Rscript.exe" "$L/dev/release/run-tests.R" "$p" "$L/$f" \
      >> "$OUT/$ct--$p--$(basename "$f" .R).txt" 2>&1 &
  done
done
wait
echo "done ob"
