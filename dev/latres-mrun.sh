#!/usr/bin/env bash
# dev/latres-m1.R and dev/latres-m4.R over their designs, N processes at a time.
source "$(dirname "$0")/latres-env.sh"
cd "$LATRES_ROOT" || exit 1
S=C:/Users/adf44/source/r/wt-latres-lib-snap
N=$1; ns1=$2; ns4=$3
for d in bern bern_long binom3 pois pois_ar9 pois_s1 pois_big bern_s05 bern_big binom5 binom10 binom20 binom10_lo pois_lo pois_vlo pois_3 nb nb_lo gamma ord4; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  ( "$RSCRIPT" dev/latres-m1.R "$S" $d 1 $ns1 dev/latres-out/m1/$d.tsv > dev/latres-out/m1/$d.log 2>&1 ) &
done
for d in gamma_log gamma_inv beta_logit weibull_log expo_log invg_log probit cloglog binom_probit pois_sqrt pois_cens gamma_big pois_ctl pois_id gamma_id gamma_invT binom_id pois_idlo; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  ( "$RSCRIPT" dev/latres-m4.R "$S" $d 1 $ns4 dev/latres-out/m4/$d.tsv > dev/latres-out/m4/$d.log 2>&1 ) &
done
wait
echo mrun-done
