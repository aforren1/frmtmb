#!/bin/bash
# predraw2's false-alarm inputs and mixture coverage on the merged build, one R process at a time.
R=/c/Users/adf44/source/r/frmtmb-wt-rel083
RUN=$R/dev/rel083-kt/run.sh
K=dev/rel083-kt/predraw2
O=dev/rel083-log/kt/predraw2
LIB=C:/Users/adf44/source/r/rellib-r22
cd $R
$RUN $R $O/kd-pois.log $K/predraw2-kdesign.R $LIB $O/predraw2-kd-pois.tsv 1:10 pois_K100_m5_mu1 pois_K100_m10_mu1.5 pois_K100_m5_mu1.5 pois_K200_m3_mu1.5 pois_K200_m5_mu1 pois_K200_m5_mu2 pois_K500_m5_mu2 pois_K50_m5_mu2; echo "kd-pois $?"
$RUN $R $O/kd-bg.log $K/predraw2-kdesign.R $LIB $O/predraw2-kd-bg.tsv 1:10 bern_K60_m20_p03 bern_K100_m10_p03 gaus_K10_m3 gaus_K50_m3 gaus_K200_m3 gaus_K500_m3; echo "kd-bg $?"
$RUN $R $O/kd-gam.log $K/predraw2-kdesign.R $LIB $O/predraw2-kd-gam.tsv 1:10 gamma_K10_m3 gamma_K50_m3 gamma_K200_m3 gamma_K500_m3; echo "kd-gam $?"
$RUN $R $O/kd-sparse.log $K/predraw2-kdesign.R $LIB $O/predraw2-kd-sparse.tsv 1:20 bern_K20_m10_p03 gamma_K20_m2 gamma_K50_m2 pois_K20_m2_mu1 pois_K50_m3_mu1; echo "kd-sparse $?"
for g in target sweep zip mix wb; do
  $RUN $R $O/meas3-$g.log $K/predraw2-measure.R $LIB $g $O/predraw2-meas3-$g-lane.rds; echo "meas3 $g $?"
done
$RUN $R $O/mc42-g.log $K/predraw2-mixcov.R $LIB mixture_gaussian:boundary 1:40 1000 $O/predraw2-mc42-g-lane.tsv; echo "mc42 g $?"
$RUN $R $O/mc42-p.log $K/predraw2-mixcov.R $LIB mixture_poisson:boundary 1:40 1000 $O/predraw2-mc42-p-lane.tsv; echo "mc42 p $?"
$RUN $R $O/mc42-c.log $K/predraw2-mixcov.R $LIB mixture_cumulative:interior 1:40 1000 $O/predraw2-mc42-c-lane.tsv; echo "mc42 c $?"
$RUN $R $O/fa-sum.txt $K/predraw2-fa-sum.R; echo "fa-sum $?"
PREDRAW2_MEAS=meas3 $RUN $R $O/meas3-cmp.txt $K/predraw2-meas-cmp.R target sweep zip mix wb; echo "meas-cmp $?"
$RUN $R $O/mc42-sum.txt $K/predraw2-mixcov-sum.R $O/predraw2-mc42; echo "mixcov-sum $?"
echo ALLDONE
