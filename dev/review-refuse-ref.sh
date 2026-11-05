#!/usr/bin/env bash
# Reference BLAS, every gate on, one process per file, final tree: the 69 files
# of the OpenBLAS list plus the gated cdf file and the plain mm and prior files.
set -u
W="$(pwd)"
OB="$W/dev/review-refuse-out/ob"
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp}"
export TREE="$W"
export FRMTMB_STAN_CACHE="$W/dev/review-refuse-out/stan-cache"
cp "$OB/list.txt" "$OB/list-ref.txt"
for f in tests/testthat/test-brms-lik-cdfs.R tests/testthat/test-prior-compat.R tests/testthat/test-prior-one-slot.R tests/testthat/test-brms-priors.R tests/testthat/test-brms-formula-priors.R tests/testthat/test-compat.R; do
  grep -q "$f" "$OB/list-ref.txt" || echo "frmtmb $f" >> "$OB/list-ref.txt"
done
for f in extensions/frmtmb.sample/tests/testthat/test-reparam.R extensions/frmtmb.sample/tests/testthat/test-simulators.R; do
  grep -q "$f" "$OB/list-ref.txt" || echo "frmtmb.sample $f" >> "$OB/list-ref.txt"
done
bash "$OB/dev/ciharden-scan.sh" reviewref ref C:/Users/adf44/source/r/wt-refuse-lib 4 "$OB/list-ref.txt" > "$OB/scan-ref.log" 2>&1
