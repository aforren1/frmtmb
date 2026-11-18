#!/usr/bin/env bash
# rel084-kt: nbrobust key tables on rellib-r23, reference BLAS only (no
# OpenBLAS R build in this tree), one R process at a time.
source /c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/env.sh
cd "$ROOT" || exit 1
K=dev/rel084-kt
O=dev/rel084-log/kt/nbrobust-out
mkdir -p $O
export OPENBLAS_NUM_THREADS=1
krun nbrobust-jac.log $K/nbrobust-jac.R $O/jac.tsv
for case in nb200 nb300; do for cse in warning ignore; do
  "$RS" $K/nbrobust-repro.R lib $case 0 20 $cse > $O/lib-ref-$case-$cse.txt 2>&1
done; done
krun nbrobust-repro-sum.log dev/nbrobust-repro-sum.R $O
krun nbrobust-start-ignore.log $K/nbrobust-start.R lib 30 ignore
krun nbrobust-start-warning.log $K/nbrobust-start.R lib 30 warning
krun nbrobust-review-fit.log $K/nbrobust-review-fit.R lib
krun nbrobust-review-fit300.log $K/nbrobust-review-fit300.R lib
echo ALLDONE > $KLOG/nbrobust-done.txt
