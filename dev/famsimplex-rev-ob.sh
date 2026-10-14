#!/usr/bin/env bash
# Reviewer: the lane's new test files under OpenBLAS 0.3.32 with LAPACK
# (a build another lane made; read only), one process per file.
#   bash dev/famsimplex-rev-ob.sh
# Output: dev/famsimplex-rev-ob/<file>.txt
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-famsimplex
RS=/c/Users/adf44/source/r/frmtmb-wt-ci701/dev/ciharden-out/Rob0.3.32-lapack/bin/Rscript.exe
OUT="$ROOT/dev/famsimplex-rev-ob"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true OPENBLAS_NUM_THREADS=4
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-famsimplex-lib
: "${TMP:?TMP is unset}"
mkdir -p "$OUT"
run() {
  p="$1"; f="$2"
  ( cd "$(dirname "$f")" && "$RS" "$ROOT/dev/release/run-tests.R" "$p" "$f" \
      > "$OUT/$(basename "$f" .R).txt" 2>&1 )
}
run frmtmb "$ROOT/tests/testthat/test-simplex-families.R" &
run frmtmb "$ROOT/tests/testthat/test-simplex-brms.R" &
run frmtmb.sample \
  "$ROOT/extensions/frmtmb.sample/tests/testthat/test-simplex-draws.R" &
wait
"$RS" -e 'cat(La_library(), extSoftVersion()["BLAS"], "\n")' \
  > "$OUT/blas.txt" 2>&1
grep -a "^RESULT\|^lib:" "$OUT"/*.txt
