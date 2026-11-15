#!/usr/bin/env bash
# Review of lane arspat2: the lane's two changed test files under R with
# OpenBLAS BLAS and LAPACK (dev/ciharden-openblas.sh <ver> lapack built
# it), one file per process, two at once. Usage:
#   bash dev/reviews/arspat2-review-openblas.sh <ver> [coretype]
set -u
VER="$1"; CT="${2:-}"
WT=C:/Users/adf44/source/r/frmtmb-wt-arspat2
RS="$WT/dev/ciharden-out/Rob$VER-lapack/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export OPENBLAS_NUM_THREADS=4 NOT_CRAN=true
if [ -n "$CT" ]; then export OPENBLAS_CORETYPE="$CT"; fi
tag="ob$VER${CT:+-$CT}"
i=0
for pf in frmtmb:tests/testthat/test-autocor-se-mi.R \
          frmtmb.spatial:extensions/frmtmb.spatial/tests/testthat/test-sar.R; do
  p="${pf%%:*}"; f="${pf#*:}"; i=$((i + 1))
  ( cd "$WT" && R_PARALLEL_PORT=$((11760 + i)) "$RS" dev/arspat2-test1.R "$p" "$f" \
      > "dev/reviews/arspat2-review-out/$tag-$(basename "$f" .R).txt" 2>&1 ) &
done
wait
grep -h "RESULT\|^lib:" "$WT"/dev/reviews/arspat2-review-out/"$tag"-*.txt
