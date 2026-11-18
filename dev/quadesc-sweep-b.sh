#!/usr/bin/env bash
# Lane quadesc: the rest of sweep "a" after its driver was stopped for
# time, as two serial chains (two R processes at most). Each chain starts
# when one of the two sweep-"a" jobs still running ends: the lane on
# OpenBLAS 0.3.26 and 0.3.32, then 0.83.0 without the Laplace reference.
set -u
cd C:/Users/adf44/source/r/frmtmb-wt-quadesc
export QUADESC_NOREF=1
run1() {
  QUADESC_REF=$1 bash dev/quadesc-run.sh "$3" default "sw-a-$2-$3-$4" \
    dev/quadesc-probe.R "$2" "$4" 0:20
}
ended() { grep -q "^END" "dev/quadesc-log/sw-a-lib-$1.txt" 2>/dev/null; }
(
  until ended 0.3.32-weak; do sleep 15; done
  run1 lap lib 0.3.26 symsmall
  run1 none base ref symsmall
  run1 none base 0.3.32 symsmall
  run1 none base 0.3.26 symsmall
  run1 none base 0.3.26 alphax
) &
(
  until ended 0.3.26-alphax; do sleep 15; done
  run1 lap lib 0.3.26 weak
  run1 none base ref weak
  run1 none base 0.3.32 weak
  run1 none base 0.3.26 weak
) &
wait
echo SWEEP DONE
