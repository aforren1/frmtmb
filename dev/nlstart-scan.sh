#!/usr/bin/env bash
# Run a per-family script over the 9 families, three R processes at once.
# Usage: bash dev/nlstart-scan.sh <arm> <outdir> [script, default
# dev/nlstart-scan.R]
arm="$1"
out="$2"
script="${3:-dev/nlstart-scan.R}"
cd /c/Users/adf44/source/r/frmtmb-wt-nlstart || exit 1
mkdir -p "$out"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NLSTART_ARM="$arm"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
fams="cumulative sratio cratio acat hurdle_cumulative bernoulli poisson gaussian gamma"
running=0
for f in $fams; do
  "$RS" "$script" "$f" > "$out/$f.txt" 2>&1 &
  running=$((running + 1))
  if [ "$running" -ge 3 ]; then
    wait -n
    running=$((running - 1))
  fi
done
wait
echo "SCAN DONE $arm $script"
