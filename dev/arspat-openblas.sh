#!/usr/bin/env bash
# Lane arspat: the lane's numerical test files under R with OpenBLAS
# BLAS and LAPACK (dev/ciharden-openblas.sh <ver> lapack builds it), three
# at a time, one file per process. Usage:
#   bash dev/arspat-openblas.sh <ver>
# Output: dev/arspat-out/ob<ver>-<file>.txt
set -u
VER="$1"
WT=C:/Users/adf44/source/r/frmtmb-wt-arspat
RS="$WT/dev/ciharden-out/Rob$VER-lapack/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export OPENBLAS_NUM_THREADS=4 NOT_CRAN=true
files="frmtmb:tests/testthat/test-autocor-se-mi.R
frmtmb:tests/testthat/test-autocor-loglik.R
frmtmb:tests/testthat/test-autocor.R
frmtmb:tests/testthat/test-sandwich.R
frmtmb:tests/testthat/test-autocor-seam.R
frmtmb.spatial:extensions/frmtmb.spatial/tests/testthat/test-sar.R"
for pf in $files; do
  p="${pf%%:*}"; f="${pf#*:}"
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 2; done
  ( cd "$WT" && "$RS" dev/arspat-test1.R "$p" "$f" \
      > "dev/arspat-out/ob$VER-$(basename "$f" .R).txt" 2>&1 ) &
done
wait
grep -h "RESULT\|^lib:" "$WT"/dev/arspat-out/ob"$VER"-*.txt
