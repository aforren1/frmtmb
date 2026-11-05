#!/usr/bin/env bash
# item 1: patched factor and control, one chain per process, 8 at once.
for arm in perm same; do
  for s in 11 12 13 14; do
    bash dev/review-gpgap-run.sh dev/review-gpgap-perm.R perm-$arm-$s.txt $arm $s 3000 12 &
  done
done
wait
