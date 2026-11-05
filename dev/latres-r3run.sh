#!/usr/bin/env bash
# The statistics of R3's condition over the three sets of designs.
source "$(dirname "$0")/latres-env.sh"
cd "$LATRES_ROOT" || exit 1
mkdir -p dev/latres-out/r3
run() {
  while [ "$(jobs -rp | wc -l)" -ge 7 ]; do sleep 1; done
  "$RSCRIPT" dev/latres-r3.R "$1" "$2" "$3" "$4" "dev/latres-out/r3/$1-$2.tsv" \
    > "dev/latres-out/r3/$1-$2.log" 2>&1 &
}
for d in bern bern_long binom3 pois pois_ar9 pois_s1 pois_big bern_s05 bern_big binom5 \
  binom10 binom20 binom10_lo pois_lo pois_vlo pois_3 nb nb_lo gamma ord4; do
  run lane $d 1 40
done
for d in zip zip_hi zip_lo hp zib cum3 cum5 cum7 nb_m1 nb_m5 nb_m20 nb_s1 beta_lo binom4 \
  ln_s2 ln_s05 sn_u1 sn_u10; do
  run m1b $d 3001 3040
done
for d in bb_T zoib_T glog_c glog_s5 glog_w blog_w ginv_b1 ginv_b2 ginv_b3 ginv_s5 gident_b \
  gid_s5 ig_b1 ig_b2 ig_s5 pid_b pid_20 pid_8 bid_b; do
  run m4b $d 2001 2040
done
wait
echo r3-done
