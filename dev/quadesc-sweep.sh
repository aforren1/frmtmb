#!/usr/bin/env bash
# Lane quadesc: the probe over arms, BLAS builds and cases, k = 0..20,
# two R processes at a time. The fit, its exact marginal log-likelihood,
# and the Laplace fit's (QUADESC_REF=lap).
#   bash dev/quadesc-sweep.sh <tag> "<arms>" "<vers>" "<cases>"
# Logs dev/quadesc-log/sw-<tag>-<arm>-<ver>-<case>.txt
set -u
TAG="$1"; ARMS="$2"; VERS="$3"; CASES="$4"
cd C:/Users/adf44/source/r/frmtmb-wt-quadesc
export QUADESC_NOREF=1 QUADESC_REF="${QUADESC_REF:-lap}"
n=0
for arm in $ARMS; do for ver in $VERS; do for c in $CASES; do
  bash dev/quadesc-run.sh "$ver" default "sw-$TAG-$arm-$ver-$c" \
    dev/quadesc-probe.R "$arm" "$c" 0:20 &
  n=$((n + 1))
  if [ "$n" -ge 2 ]; then wait -n; n=$((n - 1)); fi
done; done; done
wait
echo SWEEP DONE
