#!/usr/bin/env bash
# rel084-kt: sparsese key tables on rellib-r23, one R process at a time.
source /c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/env.sh
cd "$ROOT" || exit 1
K=dev/rel084-kt
L=C:/Users/adf44/source/r/rellib-r23
for s in 1 2 3; do krun sparsese-3370ref-s$s.log $K/sparsese-3370ref.R lane 4500 $s; done
for s in 1 2 3; do krun sparsese-fix1-s$s.log $K/sparsese-fix1.R $L $s; done
krun sparsese-fix1-s1-g0.log $K/sparsese-fix1.R $L 1 0
krun sparsese-fix1-s2-g0.log $K/sparsese-fix1.R $L 2 0
krun sparsese-cum4500.log $K/sparsese-review-cum4500.R lane default
echo ALLDONE > $KLOG/sparsese-done.txt
