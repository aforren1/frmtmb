#!/usr/bin/env bash
# rel084-kt: gpselimit key tables (the short check's single driver) on
# rellib-r23, then merge probes 8 and 4. One R process at a time.
source /c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/env.sh
cd "$ROOT" || exit 1
K=dev/rel084-kt
G=C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/reviews/gpselimit-short-out
O=C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-log/kt/gpselimit-out
krun gpselimit-short-a.log $K/gpselimit-short-a.R lane
krun gpselimit-cmp-base.log $K/gpselimit-short-cmp.R $G/a-base.rds $O/a-lane.rds
krun gpselimit-cmp-lane.log $K/gpselimit-short-cmp.R $G/a-lane.rds $O/a-lane.rds
krun probe-step8.log $K/probe-step8-drawdir-gpselimit.R
krun probe-step4.log $K/probe-step4-alapse-fitend.R
echo ALLDONE > $KLOG/gpselimit-done.txt
