#!/usr/bin/env bash
# Lane priors2, fix round: the prior and bound test files, one process
# each, OB_N at a time (default 1), under an OpenBLAS copy of R built by
# dev/ciharden-openblas.sh. Args: <ver> (lapack mode). Each log opens
# with the R home (the copy) and ends with its RESULT line. OB_FILES
# replaces the file list, OB_TAG suffixes the output directory.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-priors2
RS="$WT/dev/ciharden-out/Rob$1-lapack/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-priors2-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
OUT=$WT/dev/priors2-out/ob$1${OB_TAG:-}
mkdir -p "$OUT"
FILES="${OB_FILES:-test-prior-bounds-ends test-prior-mixture-theta test-predfix
test-priors-bounds-grcov test-mixture-starts test-prior-table-brms
test-priors-autocor-classes}"
for f in $FILES; do
  while [ "$(jobs -rp | wc -l)" -ge "${OB_N:-1}" ]; do sleep 1; done
  ( "$RS" -e "cat('blas:', R.home(), '\n')" > "$OUT/$f.txt" 2>&1
    "$RS" "$WT/dev/release/run-tests.R" frmtmb "$WT/tests/testthat/$f.R" \
      >> "$OUT/$f.txt" 2>&1 ) &
done
wait
grep -h "^blas\|^lib\|^RESULT" "$OUT"/*.txt
