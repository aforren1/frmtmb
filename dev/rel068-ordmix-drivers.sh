#!/usr/bin/env bash
# The merged-tree reruns the ordmix and fixes reviews asked for, on the
# release library: the ordmix false-alarm drivers (the round-1 table of
# dev/ordmix-p1-fa.sh, the review's margin designs, its B1 edge cases,
# its B3/B5 weight designs) and lane fixes' flat-warning table
# (dev/fixes-p2-table.R). Each rel068-<name>.R is the lane's script with
# the release library in .libPaths(); logs in dev/rel068-log/drivers/.
#
#   bash dev/rel068-ordmix-drivers.sh
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$ROOT/dev/rel068-log/drivers
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT/fa" "$OUT/margin" "$OUT/b1" "$OUT/b3"
cd "$ROOT" || exit 1
for m in cum2 probit2 cum_sr_clog mu2 hurdle_huz sr_acat_th gr_cum2 \
         mu_cum_sr_theta hurdle_prior mu_cum_links_theta mu_same_theta \
         plain_cum plain_probit; do
  ( "$RS" dev/rel068-falsealarm.R $m 20 > "$OUT/fa/$m.txt" 2>&1 ) &
done
for c in "cauchit 1 500" "cauchit 5 500" "cauchit 10 500" "cloglog 1 500" \
         "cloglog 5 500" "logit 1 5000" "logit 5 500" "logit 10 500" \
         "logit 10 2000" "probit 5 500"; do
  set -- $c
  ( "$RS" dev/rel068-margin.R $1 $2 $3 > "$OUT/margin/$1-$2-$3.txt" 2>&1 ) &
done
for k in hu1_only thres_only beta11 none nopred_mu1 disc_differs disc_same; do
  ( "$RS" dev/rel068-b1.R $k 10 > "$OUT/b1/$k.txt" 2>&1 ) &
done
for g in w_sum1 w_tenth; do
  ( "$RS" dev/rel068-b3.R $g > "$OUT/b3/$g.txt" 2>&1 ) &
done
( "$RS" dev/fixes-p2-table.R C:/Users/adf44/source/r/rellib-r6 \
    > "$OUT/fixes-p2-table.txt" 2>&1 ) &
wait
echo "DRIVERS DONE $(date +%T)"
