#!/usr/bin/env bash
# item 1: a second set of chains, patched and control; and the mode start where
# frm() has not collapsed.
for arm in perm same; do
  for s in 31 32 33 34; do
    bash dev/review-gpgap-run.sh dev/review-gpgap-perm.R perm-$arm-$s.txt $arm $s 3000 12 &
  done
done
REV_S=3e-3 R_PARALLEL_PORT=12411 bash dev/review-gpgap-run.sh dev/review-gpgap-draw3.R draw3-fit-3e-3.txt frm_fit 11 1000 10 &
REV_S=2e-3 R_PARALLEL_PORT=12477 bash dev/review-gpgap-run.sh dev/review-gpgap-draw3.R draw3-fit-2e-3.txt frm_fit 11 1000 10 &
REV_S=0.1 R_PARALLEL_PORT=12533 bash dev/review-gpgap-run.sh dev/review-gpgap-draw3.R draw3-fit-0.1.txt frm_fit 11 1000 10 &
wait
