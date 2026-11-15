#!/usr/bin/env bash
# Short check of lane priors2's fix round, item 6: the named test files on
# the lane library, one file per R process, 2 at a time (machine notice);
# then the two new files under the reviewer's OpenBLAS 0.3.32 LAPACK copy
# with OPENBLAS_CORETYPE=Haswell (no build here).
# bash dev/reviews/priors2-short-tests.sh
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-priors2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OB="$WT/dev/ciharden-out/Rob0.3.32-lapack/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-priors2-lib
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
OUT=$WT/dev/reviews/priors2-short-tests-out
mkdir -p "$OUT"
JOBS="frmtmb:tests/testthat/test-prior-bounds-ends.R
frmtmb:tests/testthat/test-prior-mixture-theta.R
frmtmb:tests/testthat/test-prior-table-brms.R
frmtmb:tests/testthat/test-mixture-starts.R
frmtmb:tests/testthat/test-conditions-census.R
frmtmb.sample:extensions/frmtmb.sample/tests/testthat/test-default-priors-autocor.R
frmtmb.sample:extensions/frmtmb.sample/tests/testthat/test-default-priors-dpar.R
frmtmb.sample:extensions/frmtmb.sample/tests/testthat/test-flat-prior.R
frmtmb.sample:extensions/frmtmb.sample/tests/testthat/test-me-prior-draws.R
frmtmb.sample:extensions/frmtmb.sample/tests/testthat/test-prior-draws.R"
port=11931
for j in $JOBS; do
  while [ "$(jobs -rp | wc -l)" -ge 2 ]; do sleep 2; done
  pkg=${j%%:*}; f=${j#*:}; b=$(basename "$f" .R)
  port=$((port + 1))
  ( R_PARALLEL_PORT=$port "$RS" "$WT/dev/release/run-tests.R" "$pkg" "$WT/$f" \
      > "$OUT/$pkg-$b.txt" 2>&1 ) &
done
wait
for f in test-prior-bounds-ends test-prior-mixture-theta; do
  while [ "$(jobs -rp | wc -l)" -ge 2 ]; do sleep 2; done
  port=$((port + 1))
  ( export OPENBLAS_NUM_THREADS=4 OPENBLAS_CORETYPE=Haswell R_PARALLEL_PORT=$port
    "$OB" -e "cat('blas:', R.home(), extSoftVersion()[['BLAS']], La_library(), '\n')" > "$OUT/ob-$f.txt" 2>&1
    "$OB" "$WT/dev/release/run-tests.R" frmtmb "$WT/tests/testthat/$f.R" \
      >> "$OUT/ob-$f.txt" 2>&1 ) &
done
wait
grep -H "^blas\|^lib\|^RESULT" "$OUT"/*.txt
