#!/usr/bin/env bash
# rel084-kt: drawdir key tables on rellib-r23, one R process at a time.
source /c/Users/adf44/source/r/frmtmb-wt-rel084/dev/rel084-kt/env.sh
cd "$ROOT" || exit 1
K=dev/rel084-kt
L=C:/Users/adf44/source/r/rellib-r23
O=dev/rel084-log/kt/drawdir-out
mkdir -p $O
export DRAWDIR_LIB=$L; krun drawdir-repro.log $K/drawdir-repro.R
krun drawdir-fires12.log $K/drawdir-fires12.R slpw:12 slpw:14 slpw:15 slpw:18 slp10:10 prsp:1 prsp:10 ordrs:9 S0.1:14 S0.1:16 S0.1:24 S0.2:8
for g in gglm gglmm wb kd zip sweep target mix40 partR; do
  krun drawdir-v-$g.log $K/drawdir-verdict.R $L $O/drawdir-v-$g.tsv $g
done
bash $K/drawdir-v-sum.sh gglm gglmm wb kd zip sweep target mix40 partR > $KLOG/drawdir-v-sum.txt 2>&1
krun drawdir-sv-main.log $K/drawdir-review-survey.R $O/sv-main.tsv bfew4 bfew6 binfew xg xp ordg ordf nbglm nbglmm dsig dphi prs prb sigid
krun drawdir-sv-b.log $K/drawdir-review-survey.R $O/sv-b.tsv xg3 nbpois prsp brs ordrs sleep slp slp10 slpw
echo ALLDONE > $KLOG/drawdir-done.txt
