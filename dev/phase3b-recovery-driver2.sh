#!/usr/bin/env bash
# Punch round 1 recovery arms, on the frozen build in phase3b-lib5.
# Six workers; the machine is shared with other lanes.
cd "$(dirname "$0")/.." || exit 1
export P3B_LIB=C:/Users/adf44/source/r/phase3b-lib7
R="C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/phase3b-log/recov3-logs
run() { "$R" dev/phase3b-eam-recovery5.R "$1" "$2" "$3" \
          >> "dev/phase3b-log/recov3-logs/$1-$2-$3.log" 2>&1; }
( run cleft 1 30;    run contfix 1 30;  run collapse 1 20 ) &
( run cleft 31 60;   run contfix 31 60; run collapse 21 40 ) &
( run contdl 1 30;   run cont 1 20;     run cens 1 20 ) &
( run contdl 31 60;  run cont 21 40;    run cens 21 40 ) &
( run cleft 61 80;   run contdl 61 80;  run cens 41 60 ) &
( run contfix 61 80; run cont 41 60;    run collapse 41 60 ) &
wait
date -Iseconds > dev/phase3b-log/recov3.done
