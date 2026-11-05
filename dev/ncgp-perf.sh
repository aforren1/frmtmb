#!/usr/bin/env bash
# dev/ncgp-perf.R through the working-set guard, base and lane interleaved, one at
# a time, so the arms see the same machine. Output dev/ncgp-out/perf.txt.
L=C:/Users/adf44/source/r/wt-ncgp-lib
O=dev/ncgp-out/perf.txt
run() { powershell -NoProfile -File dev/ncgp-guard.ps1 12000 1500 dev/ncgp-perf.R "$@" 2>&1 | grep -E 'NCGP' >> "$O"; }
: > "$O"
for rep in 1 2; do
  for n in 60 200 500 1000; do
    run base base "$n" gaussian 3
    run lane "$L" "$n" gaussian 3
  done
done
for fam in poisson bernoulli; do
  for n in 60 200; do
    run base base "$n" "$fam" 3
    run lane "$L" "$n" "$fam" 3
  done
done
echo PERF-DONE >> "$O"
