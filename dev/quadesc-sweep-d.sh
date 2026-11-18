#!/usr/bin/env bash
# Lane quadesc, fix round 1: the fix-round build (a restart counts as a
# gain only when it left the point; an abandoned point carries its
# counts) over the grid of dev/quadesc-sweep-c.sh: k = 0..20, three
# cases, three BLAS builds, two R processes at a time, no Laplace
# reference. Logs sw-d-lib-<ver>-<case>.txt
set -u
cd C:/Users/adf44/source/r/frmtmb-wt-quadesc
export QUADESC_NOREF=1 QUADESC_REF=none
run1() {
  bash dev/quadesc-run.sh "$1" default "sw-d-lib-$1-$2" \
    dev/quadesc-probe.R lib "$2" 0:20
}
(
  run1 ref alphax; run1 0.3.32 alphax; run1 0.3.26 alphax; run1 ref weak
) &
(
  run1 ref symsmall; run1 0.3.32 symsmall; run1 0.3.26 symsmall
  run1 0.3.32 weak; run1 0.3.26 weak
) &
wait
echo SWEEP DONE
