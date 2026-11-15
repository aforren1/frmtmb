#!/usr/bin/env bash
# Short check of lane gpfield: test files under OpenBLAS 0.3.32 (lapack mode,
# Haswell kernels), on an R copy built in the session scratch directory.
# Args: <test file stems...>; 1 at a time.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-gpfield
RS=/c/Users/adf44/AppData/Local/Temp/2/claude/ee0bfefdf562474c/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/obroot/dev/ciharden-out/Rob0.3.32-lapack/bin/Rscript.exe
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4 OPENBLAS_CORETYPE=Haswell
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-gpfield-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
OUT=$WT/dev/reviews/gpfield-short-out/ob0332
mkdir -p "$OUT"
for f in "$@"; do
  while [ "$(jobs -rp | wc -l)" -ge 1 ]; do sleep 1; done
  ( "$RS" -e "cat('blas:', R.home(), extSoftVersion()[['BLAS']], La_library(), '\n')" > "$OUT/$f.txt" 2>&1
    "$RS" "$WT/dev/release/run-tests.R" frmtmb "$WT/tests/testthat/$f.R" >> "$OUT/$f.txt" 2>&1 ) &
done
wait
grep -h "^blas\|^lib\|^RESULT" "$OUT"/*.txt
