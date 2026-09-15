#!/bin/sh
# Run the reviewer's probe: one fresh R process per (build, mode, arm),
# three at a time.
export PATH="/c/Program Files/R/R-4.6.1/bin/x64:$PATH"
cd /c/Users/adf44/source/r/frmtmb-wt-samplegen
OUT=dev/sgrev-out
mkdir -p $OUT/probe $OUT/log
n=0
for build in BASE FIX; do
  for mode in S T U N P Q G D R RN RN2 BR CA CA2; do
    for arm in test control; do
      f="$OUT/probe/$build-$mode-$arm.rds"
      Rscript --vanilla dev/sgrev-probe.R "$build" "$mode" "$arm" "$f" \
        > "$OUT/log/$build-$mode-$arm.log" 2>&1 &
      n=$((n+1))
      if [ $((n % 3)) -eq 0 ]; then wait; fi
    done
  done
done
wait
echo "done"
ls $OUT/probe | wc -l
