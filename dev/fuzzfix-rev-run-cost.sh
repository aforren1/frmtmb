#!/usr/bin/env bash
# Reviewer: run cost configurations one process each, memory capped.
#   bash dev/fuzzfix-rev-run-cost.sh <lib|base> <label> <capMB> \
#        "<st> <ng> <m> [bypass]" ...
WT=C:/Users/adf44/source/r/frmtmb-wt-fuzzfix
lib=$1; lbl=$2; cap=$3; shift 3
cd /c/Users/adf44/source/r/frmtmb-wt-fuzzfix
for cfg in "$@"; do
  set -- $cfg
  tag=$(echo "$cfg" | tr ' ' '-')
  out="$WT/dev/fuzzfix-rev-out/cost-$tag-$lbl.txt"
  powershell -NoProfile -ExecutionPolicy Bypass -File \
    "$WT/dev/fuzzfix-rev-capped.ps1" "$cap" "$out" \
    "$WT/dev/fuzzfix-rev-quad-cost.R" "$lib" $cfg
done
