#!/usr/bin/env bash
# item 1: no harm at ordinary noise? n 40, noise 0.3, patched and control; and the
# patched factor at low noise, n 40, short.
export REV_S=0.3
for arm in perm same; do
  for s in 21 22 23 24; do
    bash dev/review-gpgap-run.sh dev/review-gpgap-perm2.R perm2-$arm-$s.txt $arm $s 1000 10 40 &
  done
done
wait
