#!/usr/bin/env bash
# Reviewer, fix round 1: run arms of a script one capped process each.
#   bash dev/fuzzfix-rev-r1-run.sh <script> <lib|base> <label> <capMB> <arm...>
WT=C:/Users/adf44/source/r/frmtmb-wt-fuzzfix
sc=$1; lib=$2; lbl=$3; cap=$4; shift 4
cd /c/Users/adf44/source/r/frmtmb-wt-fuzzfix
stem=$(basename "$sc" .R)
for arm in "$@"; do
  out="$WT/dev/fuzzfix-rev-out/${stem#fuzzfix-rev-}-$arm-$lbl.txt"
  powershell -NoProfile -ExecutionPolicy Bypass -File \
    "$WT/dev/fuzzfix-rev-capped.ps1" "$cap" "$out" "$WT/dev/$sc" "$lib" $arm
done
