#!/usr/bin/env bash
# every mutant, alone and with each extension set the change was for
cd "$(dirname "$0")/.."
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
out=dev/testiso-log/mutants-$1; rm -rf $out; mkdir -p $out
i=0
for m in none none_reg contrib_dead core_unresolved core_special core_method unresolved_drop commit_early aterm_nofeature; do
  for pre in "" "frmtmb.sample" "frmtmb.eam" "frmtmb.ode" "frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline"; do
    i=$((i+1))
    while [ "$(jobs -rp | wc -l)" -ge 6 ]; do sleep 1; done
    "$RS" dev/testiso-mutants.R $m $pre > $out/$i.txt 2>&1 &
  done
done
wait
for k in $(seq 1 $i); do grep -a -E "^mutant|^  FAILS" $out/$k.txt; done
