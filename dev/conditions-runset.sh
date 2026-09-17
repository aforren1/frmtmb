#!/usr/bin/env bash
# Lane wt-conditions: run test files one per R process, P at a time,
# through the sweep runner, on one library arm.
#   dev/conditions-runset.sh <lane|base> <P> <outdir> file...
arm=$1; P=$2; out=$3; shift 3
mkdir -p "$out"
date +%H:%M:%S > "$out/.start"
printf '%s\n' "$@" | xargs -P "$P" -I{} bash -c '
  f="{}"; log="'"$out"'/$(echo "$f" | tr "/" "_").log"
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" \
    dev/conditions-sweep-onefile.R '"$arm"' "$f" '"$out"' > "$log" 2>&1
  l=$(grep -h "^BLOCKS" "$log")
  if [ -z "$l" ]; then l="NO RESULT LINE  $f"; fi
  echo "$l"
' > "$out/summary.txt"
date +%H:%M:%S > "$out/.end"
