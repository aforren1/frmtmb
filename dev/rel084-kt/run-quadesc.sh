#!/usr/bin/env bash
# rel084-kt: quadesc key tables on rellib-r23, reference BLAS only (no
# OpenBLAS R build in the rel084 tree), one R process at a time. The lane's
# dev/quadesc-run.sh set OPENBLAS_NUM_THREADS=4 and unset OPENBLAS_CORETYPE;
# kept here (no effect on reference BLAS).
source /c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/env.sh
cd "$ROOT" || exit 1
K=dev/rel084-kt
Q=dev/rel084-log/kt/quadesc-log
export OPENBLAS_NUM_THREADS=4; unset OPENBLAS_CORETYPE
for c in alphax symsmall weak; do
  QUADESC_NOREF=1 QUADESC_REF=none "$RS" $K/quadesc-probe.R lib $c 0:20 > $Q/sw-d-lib-ref-$c.txt 2>&1
  echo "exit $?" >> $Q/sw-d-lib-ref-$c.txt
done
krun quadesc-sweep-sum.log $K/quadesc-sweep-sum.R d
"$RS" $K/quadesc-veto-ulp.R lib > $Q/veto-lib-ref.txt 2>&1
echo "exit $?" >> $Q/veto-lib-ref.txt
echo ALLDONE > $KLOG/quadesc-done.txt
