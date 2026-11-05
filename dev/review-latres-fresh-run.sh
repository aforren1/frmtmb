#!/bin/bash
cd /c/Users/adf44/source/r/frmtmb-wt-latres
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for d in f_pois05 f_pois_ar8 f_pois_short f_pois_long f_bin2 f_bin6 f_bin30 f_probit f_cloglog f_geom f_nb1 f_nb8 f_gam1 f_gam2s1 f_gam20 f_wei2 f_exp f_beta2 f_beta20 f_zinb f_hnb f_cumpro f_ln1; do
  echo "$d 7001 7030 dev/review-latres-out/fresh/$d.tsv"
done | xargs -P 7 -L 1 "$R" dev/review-latres-fresh.R > dev/review-latres-out/fresh.log 2>&1
for d in f_pois05 f_pois_ar8 f_pois_short f_pois_long f_bin2 f_bin6 f_bin30 f_probit f_cloglog f_geom f_nb1 f_nb8 f_gam1 f_gam2s1 f_gam20 f_wei2 f_exp f_beta2 f_beta20 f_zinb f_hnb f_cumpro f_ln1; do
  echo "fresh $d 7001 7030 dev/review-latres-out/tk/fresh-$d.tsv"
done > dev/review-latres-out/tk.jobs
for d in bern bern_long binom3 pois pois_ar9 pois_s1 pois_big bern_s05 bern_big binom5 binom10 binom20 binom10_lo pois_lo pois_vlo pois_3 nb nb_lo gamma ord4; do
  echo "lane $d 1 40 dev/review-latres-out/tk/lane-$d.tsv"
done >> dev/review-latres-out/tk.jobs
for d in zip zip_hi zip_lo hp zib cum3 cum5 cum7 nb_m1 nb_m5 nb_m20 nb_s1 beta_lo binom4 ln_s2 ln_s05 sn_u1 sn_u10; do
  echo "m1b $d 3001 3040 dev/review-latres-out/tk/m1b-$d.tsv"
done >> dev/review-latres-out/tk.jobs
for d in ginv_b1 ginv_b2 ginv_b3 ginv_s5 gident_b gid_s5 ig_b1 ig_b2 ig_s5 pid_b pid_20 pid_8 bid_b glog_c glog_s5 glog_w blog_w gid_l bid_l ginv_l pid_l bb_T zoib_T; do
  echo "m4b $d 2001 2040 dev/review-latres-out/tk/m4b-$d.tsv"
done >> dev/review-latres-out/tk.jobs
xargs -P 7 -L 1 "$R" dev/review-latres-tk.R < dev/review-latres-out/tk.jobs > dev/review-latres-out/tk.log 2>&1
