#!/usr/bin/env bash
# Punch round 1, B1: the reviewer's dev/ordmix-rev-falsealarm.R, every
# model it defines (its own driver leaves out two), one process each,
# on the lane build. Logs dev/ordmix-p1-fa-log/<model>.txt; the
# reviewer's logs in dev/ordmix-rev-fa-log/ stay as the before record.
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordmix"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p "$ROOT/dev/ordmix-p1-fa-log"
for m in cum2 probit2 cum_sr_clog mu2 hurdle_huz sr_acat_th gr_cum2 \
         mu_cum_sr_theta hurdle_prior mu_cum_links_theta mu_same_theta \
         plain_cum plain_probit; do
  ( "$RS" "$ROOT/dev/ordmix-rev-falsealarm.R" $m 20 \
      > "$ROOT/dev/ordmix-p1-fa-log/$m.txt" 2>&1 ) &
done
wait
