#!/usr/bin/env bash
# Reviewer: run the quadrature shape arms one process each, memory capped.
#   bash dev/fuzzfix-rev-run-shapes.sh <lib|base> <label> <capMB> <arm...>
WT=C:/Users/adf44/source/r/frmtmb-wt-fuzzfix
lib=$1; lbl=$2; cap=$3; shift 3
cd /c/Users/adf44/source/r/frmtmb-wt-fuzzfix
for arm in "$@"; do
  out="$WT/dev/fuzzfix-rev-out/shape-$arm-$lbl${NG:+-ng$NG}.txt"
  powershell -NoProfile -ExecutionPolicy Bypass -File \
    "$WT/dev/fuzzfix-rev-capped.ps1" "$cap" "$out" \
    "$WT/dev/fuzzfix-rev-quad-shapes.R" "$lib" "$arm" ${NG:-8}
done
