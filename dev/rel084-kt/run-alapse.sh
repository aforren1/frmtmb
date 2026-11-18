#!/usr/bin/env bash
# rel084-kt: alapse key tables on rellib-r23, one R process at a time.
source /c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/env.sh
K=dev/rel084-kt
krun alapse-rq.log $K/alapse-rq.R lane 1:25
krun alapse-rqsum.log dev/alapse-rqsum2.R dev/rel084-log/kt/alapse-rq.log
krun alapse-qglob.log $K/alapse-review-qglob.R lane "c(2,3,4,12,20,22,101:140)"
krun alapse-search.log $K/alapse-short-search.R "c(109, 22, 140, 2, 3, 4, 12, 20, 101, 102, 601:610)"
krun alapse-multistart.log $K/alapse-review-multistart.R lane 1:20
krun alapse-fits.log $K/alapse-short-fits.R
krun alapse-bounds.log $K/alapse-short-bounds.R lane
krun alapse-over.log $K/alapse-short-over.R
krun alapse-review-rq.log $K/alapse-review-rq.R lane
krun alapse-ties.log $K/alapse-short-ties.R
krun alapse-sync.log $K/alapse-short-sync.R
echo ALLDONE > $KLOG/alapse-done.txt
