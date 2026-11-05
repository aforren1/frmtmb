#!/usr/bin/env bash
# The seven extensions no lane changed, every test file, on the 0.78.0
# core (r079-floor-others: rellib-r16 alone) against the merged core
# (r079-ref: rellib-r17), every gate on: the RESULT line of each file,
# and each package's totals.
#   bash dev/rel079-others-cmp.sh > dev/rel079-log/others-cmp.txt
cd /c/Users/adf44/source/r/frmtmb-wt-rel079 || exit 1
L=dev/ciharden-log
n=0; d=0
while read -r p f; do
  b=$(basename "$f" .R)
  o=$(grep -h "^RESULT" "$L/r079-floor-others-suite/$p--$b.txt" | tail -1)
  m=$(grep -h "^RESULT" "$L/r079-ref-suite/$p--$b.txt" | tail -1)
  n=$((n + 1))
  [ "$o" = "$m" ] || { d=$((d + 1)); echo "DIFFERS $p $b: 0.78.0 core: $o | merged: $m"; }
done < dev/rel079-log/floor-others.jobs
echo "lib: lines of the old-core run naming frmtmb 0.78.0 from rellib-r16:" \
  "$(grep -h '^lib:' $L/r079-floor-others-suite/*.txt | grep -c 'frmtmb 0.78.0 C:/Users/adf44/source/r/rellib-r16')" \
  "of $(ls $L/r079-floor-others-suite/*.txt | wc -l)"
for p in frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.spatial frmtmb.survival; do
  for r in floor-others ref; do
    s=$(grep -h "^RESULT" $L/r079-$r-suite/$p--*.txt | awk '{for(i=1;i<=NF;i++){split($i,a,"="); if(a[1]=="pass")p+=a[2]; if(a[1]=="fail")f+=a[2]; if(a[1]=="err")e+=a[2]; if(a[1]=="skip")s+=a[2]}; n++} END {print n " files, pass " p ", fail " f ", err " e ", skip " s}')
    printf '%s %s: %s\n' "$p" "$r" "$s"
  done
done
echo "files compared: $n, RESULT lines that differ: $d"
