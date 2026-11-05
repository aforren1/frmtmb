#!/usr/bin/env bash
# item 1: the five arms at two Stan seeds, two arms at once, ports apart.
set -u
p=11930
for s in 11 22; do
  for r in brms brms_sort frm frm_warm frm_fit; do
    while [ "$(jobs -rp | wc -l)" -ge 2 ]; do sleep 2; done
    p=$((p + 37))
    R_PARALLEL_PORT=$p bash dev/review-gpgap-run.sh dev/review-gpgap-draw.R draw-$r-$s.txt $r $s 3000 12 &
  done
done
wait
