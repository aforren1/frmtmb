#!/usr/bin/env bash
# The seven extensions no lane changed, every test file, on the 0.77.0
# core (r078-floor-others: rellib-r15 alone) against the merged core
# (r078-ref: rellib-r16), every gate on: the RESULT line of each file,
# and each package's totals.
#   bash dev/rel078-others-cmp.sh > dev/rel078-log/others-cmp.txt
cd /c/Users/adf44/source/r/frmtmb-wt-rel078 || exit 1
L=dev/ciharden-log
n=0; d=0
while read -r p f; do
  b=$(basename "$f" .R)
  o=$(grep -h "^RESULT" "$L/r078-floor-others-suite/$p--$b.txt" | tail -1)
  m=$(grep -h "^RESULT" "$L/r078-ref-suite/$p--$b.txt" | tail -1)
  n=$((n + 1))
  [ "$o" = "$m" ] || { d=$((d + 1)); echo "DIFFERS $p $b: 0.77.0 core: $o | merged: $m"; }
done < dev/rel078-log/floor-others.jobs
echo "lib: lines of the old-core run naming frmtmb 0.77.0 from rellib-r15:" \
  "$(grep -h '^lib:' $L/r078-floor-others-suite/*.txt | grep -c 'frmtmb 0.77.0 C:/Users/adf44/source/r/rellib-r15')" \
  "of $(ls $L/r078-floor-others-suite/*.txt | wc -l)"
for p in frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.spatial frmtmb.survival; do
  for r in floor-others ref; do
    s=$(grep -h "^RESULT" $L/r078-$r-suite/$p--*.txt | awk '{for(i=1;i<=NF;i++){split($i,a,"="); if(a[1]=="pass")p+=a[2]; if(a[1]=="fail")f+=a[2]; if(a[1]=="err")e+=a[2]; if(a[1]=="skip")s+=a[2]}; n++} END {print n " files, pass " p ", fail " f ", err " e ", skip " s}')
    printf '%s %s: %s\n' "$p" "$r" "$s"
  done
done
echo "files compared: $n, RESULT lines that differ: $d"
