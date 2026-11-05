#!/usr/bin/env bash
# The final pass of fix round 1: mutants, OpenBLAS, then one check per package.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for m in control class_drop no_warning thr_half thr_double thr_x5 thr_x30 thr_fixed q_max q_med nug_x10 sd1 no_trycatch no_prior first_only one_block_loop msg_no_moves msg_k_promise; do
  while [ "$(jobs -r | wc -l)" -ge 4 ]; do sleep 2; done
  ( "$RS" dev/gpgap-mut.R "$m" > "dev/gpgap-out/mut/$m.txt" 2>&1 ) &
done
wait
rm -rf tests/testthat/_problems
bash dev/gpgap-ob.sh
rm -rf dev/ciharden-out extensions/frmtmb.sample/tests/testthat/_problems tests/testthat/_problems
for p in . extensions/frmtmb.spatial extensions/frmtmb.sample; do
  bash dev/gpgap-check.sh "$p" > "dev/gpgap-out/checkrun-$(basename "$p").txt" 2>&1
done
