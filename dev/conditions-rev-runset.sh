#!/usr/bin/env bash
# Reviewer, lane wt-conditions: every test file of all 8 packages, one
# per R process, P at a time, on one library arm.
#   dev/conditions-rev-runset.sh <lane|base> <P> <outdir>
arm=$1; P=$2; out=$3
mkdir -p "$out"
date +%H:%M:%S > "$out/.start"
ls tests/testthat/test-*.R extensions/*/tests/testthat/test-*.R |
xargs -P "$P" -I{} bash -c '
  f="{}"; log="'"$out"'/$(echo "$f" | tr "/" "_").log"
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" \
    dev/conditions-rev-sweep-onefile.R '"$arm"' "$f" '"$out"' > "$log" 2>&1
  l=$(grep -h "^BLOCKS" "$log")
  if [ -z "$l" ]; then l="NO RESULT LINE  $f"; fi
  echo "$l"
' > "$out/summary.txt"
date +%H:%M:%S > "$out/.end"
