#!/usr/bin/env bash
# Final runs, second track: OpenBLAS 0.3.26 and 0.3.32 with LAPACK,
# the gated files that use importance =, and p124 and p069 against the baseline.
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp || exit 1
L=dev/impsamp-out
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LIB=C:/Users/adf44/source/r/wt-impsamp-lib
{ for f in importance grad-verdict se-exact-hessian rl-example covstruct-seam; do
    echo "frmtmb tests/testthat/test-$f.R"; done
  for f in factorization rlddm-ndt session surface; do
    echo "frmtmb.learn extensions/frmtmb.learn/tests/testthat/test-$f.R"; done
  echo "frmtmb.spatial extensions/frmtmb.spatial/tests/testthat/test-car.R"
} > $L/ob-list.txt
NOT_CRAN=true bash dev/ciharden-scan.sh impsamp-ob26 ob0.3.26-lapack-t4 "$LIB" 2 $L/ob-list.txt > $L/ob26-scan.txt 2>&1
NOT_CRAN=true bash dev/ciharden-scan.sh impsamp-ob32 ob0.3.32-lapack-t4 "$LIB" 2 $L/ob-list.txt > $L/ob32-scan.txt 2>&1
echo ob-done > $L/final-b-ob.txt
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE="$(cd dev/stan-cache && pwd -W)"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
mkdir -p $L/gated
gated() {
  "$R" dev/impsamp-testfile.R "$LIB" "$1" "$2" > "$L/gated/$(basename "$2").txt" 2>&1
}
gated frmtmb tests/testthat/test-me.R &
gated frmtmb tests/testthat/test-rl-example.R &
wait
gated frmtmb.spatial extensions/frmtmb.spatial/tests/testthat/test-brms-likelihood-car.R &
gated frmtmb.sample extensions/frmtmb.sample/tests/testthat/test-loo.R &
wait
gated frmtmb tests/testthat/test-importance.R &
gated frmtmb tests/testthat/test-bcm-simple.R &
wait
unset FRMTMB_BRMS_FIT_TESTS
echo gated-done > $L/final-b-gated.txt
rm -f dev/impsamp-equiv-out/p124.tsv dev/impsamp-equiv-out/p069.tsv
FRMTMB_EQUIV_LIB="$LIB;C:/Users/adf44/source/r/rellib-r16" bash dev/equiv-run.sh 2 dev/impsamp-equiv-out > $L/equiv-run.txt 2>&1
echo equiv-done > $L/final-b-equiv.txt
