#!/usr/bin/env bash
# Punch round 2b, on the final lane build, in turn: the core suite (one
# file per process, every gate on), then the review's two drivers
# (dev/ordmix-rev3-b3.R w_sum1 and w_tenth, dev/ordmix-rev3-wdegen.R)
# and the criterion measurements (dev/ordmix-p2-crit.sh, the margin,
# 160 and 140 populations).
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordmix"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export TMP="C:/Users/adf44/AppData/Local/Temp/1" TEMP="C:/Users/adf44/AppData/Local/Temp/1"
bash "$ROOT/dev/ordmix-run-par.sh" lane p2bcore "$ROOT/dev/ordmix-p2b-jobs-core.txt" 20
mkdir -p "$ROOT/dev/ordmix-p2b-log"
( "$RS" "$ROOT/dev/ordmix-rev3-b3.R" w_sum1 > "$ROOT/dev/ordmix-p2b-log/b3-w_sum1.txt" 2>&1 ) &
( "$RS" "$ROOT/dev/ordmix-rev3-b3.R" w_tenth > "$ROOT/dev/ordmix-p2b-log/b3-w_tenth.txt" 2>&1 ) &
( "$RS" "$ROOT/dev/ordmix-rev3-wdegen.R" > "$ROOT/dev/ordmix-p2b-log/wdegen.txt" 2>&1 ) &
wait
bash "$ROOT/dev/ordmix-p2-crit.sh"
echo FINAL DONE
