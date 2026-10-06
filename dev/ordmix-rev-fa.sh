#!/usr/bin/env bash
# Reviewer of lane ordmix: every model of dev/ordmix-rev-falsealarm.R,
# one process each. Logs dev/ordmix-rev-fa-log/<model>.txt
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordmix"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for m in cum2 probit2 cum_sr_clog mu2 hurdle_huz sr_acat_th gr_cum2 mu_cum_sr_theta hurdle_prior plain_cum plain_probit; do
  ( "$RS" "$ROOT/dev/ordmix-rev-falsealarm.R" $m 20 > "$ROOT/dev/ordmix-rev-fa-log/$m.txt" 2>&1 ) &
done
wait
