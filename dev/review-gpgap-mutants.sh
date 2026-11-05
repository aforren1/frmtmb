#!/usr/bin/env bash
# Runs the instrument's mutants, 5 at once.
set -u
for id in p094 p022; do
for m in none kern_1.1 kern_1.01 kern_1.001 kern_rq nug_1e-5 nug_2e-6 nug_1.1e-6 nug_1e-7 tr_lscale_1.01 tr_lscale_1.0001 tr_sdgp_1.01 tr_sdgp_1.0001 tr_noorder; do
  while [ "$(jobs -rp | wc -l)" -ge 5 ]; do sleep 1; done
  bash dev/review-gpgap-run.sh dev/review-gpgap-mutant.R mut-$id-$m.txt $id $m &
done
done
wait
