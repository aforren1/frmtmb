#!/usr/bin/env bash
# Lane quadesc: the lane's second build (a restart the optimizer abandoned
# is kept only when no restart that ended on its own gained), k = 0..20,
# three cases, three BLAS builds, two R processes at a time, no Laplace
# reference (its exact value does not move with k:
# dev/quadesc-log/sw-a-lib-*.txt). Logs sw-c-lib-<ver>-<case>.txt
set -u
cd C:/Users/adf44/source/r/frmtmb-wt-quadesc
export QUADESC_NOREF=1 QUADESC_REF=none
run1() {
  bash dev/quadesc-run.sh "$1" default "sw-c-lib-$1-$2" \
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
