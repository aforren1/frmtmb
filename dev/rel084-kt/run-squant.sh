#!/usr/bin/env bash
# rel084-kt: squant key tables on rellib-r23, one R process at a time.
# Not run: dev/reviews/squant-short-compat.R and squant-short-gp.R (gp()
# fits at 150 rows over 19 minutes each), the SE study, compat counts.
source /c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/env.sh
cd "$ROOT" || exit 1
K=dev/rel084-kt
krun squant-qgam.log $K/squant-qgam.R
krun squant-ulps.log $K/squant-ulps.R
krun squant-tapegrad.log $K/squant-tapegrad.R C:/Users/adf44/source/r/rellib-r23
krun squant-tapegrad2.log $K/squant-tapegrad2.R
krun squant-short-tape.log $K/squant-short-tape.R
krun squant-review-tape4.log $K/squant-review-tape4.R
krun squant-review-fits.log $K/squant-review-fits.R
echo ALLDONE > $KLOG/squant-done.txt
