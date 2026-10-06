#!/usr/bin/env bash
# Punch round 2, B3: every config of dev/ordmix-p2-crit.R, one process
# each, all at once. Logs dev/ordmix-p2-crit-log/<config>.txt
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordmix"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p "$ROOT/dev/ordmix-p2-crit-log"
for c in margin_logit_1_5000 margin_logit_5_500 margin_probit_5_500 \
         margin_cloglog_5_500 margin_logit_10_500 margin_logit_10_2000 \
         margin_cauchit_10_500 margin_cauchit_1_500 margin_cauchit_5_500 \
         margin_cloglog_1_500 margin_probit_10_500 margin_cloglog_10_500 \
         rev_cum_cum_300 rev_cum_cum_500 rev_cum_sr_300 rev_cum_sr_500 \
         id_cum2 id_probit2 id_cum_sr_clog id_mu2 id_hurdle_huz \
         id_sr_acat_th id_gr_cum2; do
  ( "$RS" "$ROOT/dev/ordmix-p2-crit.R" $c > "$ROOT/dev/ordmix-p2-crit-log/$c.txt" 2>&1 ) &
done
wait
