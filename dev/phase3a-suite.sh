#!/usr/bin/env bash
# One test file per R process, for one package source tree and one arm.
# Usage: phase3a-suite.sh <pkg source dir> <arm: base|lane> <label>
# Writes dev/phase3a-log/suite-<label>.tsv (one row per test_that block)
# and suite-<label>/<file>.txt, and prints one FILE line per file.
set -u
src="$1"; arm="$2"; label="$3"
W=C:/Users/adf44/source/r/frmtmb-wt-phase3a
out="$W/dev/phase3a-log/suite-$label"
tsv="$W/dev/phase3a-log/suite-$label.tsv"
rm -rf "$out" "$tsv"; mkdir -p "$out"
ran=0
for f in "$src"/tests/testthat/test-*.R; do
  b=$(basename "$f")
  "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" "$W/dev/phase3a-runtest.R" \
    "$src" "$b" "$arm" "$tsv" > "$out/$b.txt" 2>&1
  line=$(grep "^FILE" "$out/$b.txt")
  if [ -z "$line" ]; then line="FILE $b ARM $arm DIED"; fi
  echo "$line"
  ran=$((ran + 1))
done
echo "SUITE $label files $ran"
