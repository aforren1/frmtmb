#!/usr/bin/env bash
# rel084-kt: slopetrig key tables on rellib-r23, one R process at a time.
source /c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/env.sh
cd "$ROOT" || exit 1
K=dev/rel084-kt
O=dev/rel084-log/kt/slopetrig-out
mkdir -p $O
krun slopetrig-repro.log $K/slopetrig-repro.R lane
krun slopetrig-lme4.log $K/slopetrig-lme4.R lane
krun slopetrig-grby.log $K/slopetrig-grby.R
krun slopetrig-probe.log $K/slopetrig-review-probe.R
krun slopetrig-cover.log $K/slopetrig-cover.R
krun slopetrig-calib-gau.log $K/slopetrig-calib.R gau lane $O/calib-gau-lane.tsv
krun slopetrig-calib-pois.log $K/slopetrig-calib.R pois lane $O/calib-pois-lane.tsv
krun slopetrig-calib2-gau.log $K/slopetrig-calib2.R gau lane $O/calib2-gau-lane.tsv
krun slopetrig-calib2-pois.log $K/slopetrig-calib2.R pois lane $O/calib2-pois-lane.tsv
krun slopetrig-calib-sum.log dev/slopetrig-calib-sum.R $O/calib-gau-lane.tsv $O/calib-pois-lane.tsv
krun slopetrig-calib2-sum.log dev/slopetrig-calib-sum.R $O/calib2-gau-lane.tsv $O/calib2-pois-lane.tsv
krun slopetrig-cor.log $K/slopetrig-cor.R lane
krun slopetrig-guard.log $K/slopetrig-guard.R lane
echo ALLDONE > $KLOG/slopetrig-done.txt
