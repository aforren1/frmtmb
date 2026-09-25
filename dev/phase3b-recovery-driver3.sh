#!/usr/bin/env bash
# Punch round 1, second launch, on the frozen final build in
# phase3b-lib7: the arms the first launch lost when its cleft workers
# were stopped (every cleft fit on phase3b-lib5 ended in a NaN gradient,
# fixed since), and cleft and cens on the final build.
cd "$(dirname "$0")/.." || exit 1
export P3B_LIB=C:/Users/adf44/source/r/phase3b-lib7
R="C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/phase3b-log/recov4-logs
run() { "$R" dev/phase3b-eam-recovery6.R "$1" "$2" "$3" \
          >> "dev/phase3b-log/recov4-logs/$1-$2-$3.log" 2>&1; }
( run cleft 1 27;    run contfix 1 30 ) &
( run cleft 28 54;   run contfix 31 60 ) &
( run cleft 55 80;   run collapse 1 20 ) &
( run cens 1 30;     run collapse 21 40;  run contdl 61 80 ) &
( run cens 31 60;    run collapse 41 60 ) &
wait
date -Iseconds > dev/phase3b-log/recov4.done
