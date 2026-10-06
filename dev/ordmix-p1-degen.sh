#!/usr/bin/env bash
# Every config of dev/ordmix-p1-degen.R, one process each, all at once.
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordmix"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p "$ROOT/dev/ordmix-p1-degen-log"
for c in rev_cum_cum_300 rev_cum_cum_500 rev_cum_sr_300 rev_cum_sr_500 \
         id_cum2 id_probit2 id_cum_sr_clog id_mu2 id_hurdle_huz \
         id_sr_acat_th id_gr_cum2; do
  ( "$RS" "$ROOT/dev/ordmix-p1-degen.R" $c > "$ROOT/dev/ordmix-p1-degen-log/$c.txt" 2>&1 ) &
done
wait
