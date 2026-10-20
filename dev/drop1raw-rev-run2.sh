#!/usr/bin/env bash
# Reviewer batch 2: probes, step/drop1 cases, MAP, memory (both arms).
ROOT=/c/Users/adf44/source/r/frmtmb-wt-drop1raw
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/drop1raw-rev-out"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export NOT_CRAN=true
cd "$ROOT"
"$RS" dev/drop1raw-rev-probe.R lane > "$OUT/probe-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-probe.R base > "$OUT/probe-base.txt" 2>&1 &
"$RS" dev/drop1raw-rev-step.R lane > "$OUT/step-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-map.R lane > "$OUT/map-lane.txt" 2>&1 &
"$RS" dev/drop1raw-memory.R lane > "$OUT/memory-lane.txt" 2>&1 &
"$RS" dev/drop1raw-memory.R base > "$OUT/memory-base.txt" 2>&1 &
wait
echo done
