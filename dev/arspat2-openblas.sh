#!/usr/bin/env bash
# Lane arspat2: the lane's numerical test files under R with OpenBLAS
# BLAS and LAPACK (dev/ciharden-openblas.sh <ver> lapack builds it), three
# at a time, one file per process. Usage:
#   bash dev/arspat2-openblas.sh <ver> [tag]
# OPENBLAS_CORETYPE from the caller's environment reaches the R processes
# (Haswell, as the reviewer ran it). Output:
# dev/arspat2-out/ob<tag>-<file>.txt, tag defaulting to <ver>.
set -u
VER="$1"
TAG="${2:-$VER}"
WT=C:/Users/adf44/source/r/frmtmb-wt-arspat2
RS="$WT/dev/ciharden-out/Rob$VER-lapack/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export OPENBLAS_NUM_THREADS=4 NOT_CRAN=true
files="frmtmb:tests/testthat/test-autocor-se-mi.R
frmtmb:tests/testthat/test-autocor-seam.R
frmtmb:tests/testthat/test-autocor.R
frmtmb:tests/testthat/test-mi-families.R
frmtmb.spatial:extensions/frmtmb.spatial/tests/testthat/test-sar.R
frmtmb.spatial:extensions/frmtmb.spatial/tests/testthat/test-car.R"
for pf in $files; do
  p="${pf%%:*}"; f="${pf#*:}"
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 2; done
  ( cd "$WT" && "$RS" dev/arspat2-test1.R "$p" "$f" \
      > "dev/arspat2-out/ob$TAG-$(basename "$f" .R).txt" 2>&1 ) &
done
wait
grep -h "RESULT\|^lib:" "$WT"/dev/arspat2-out/ob"$TAG"-*.txt
