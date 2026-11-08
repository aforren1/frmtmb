#!/usr/bin/env bash
# Run dev/superassign-trace.R on each "<pkg> <test file>" line of a job
# list, base arm then lane arm, ONE R process at a time.
#
#   bash dev/superassign-trace-run.sh <job list>
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/superassign-trace-out"
mkdir -p "$OUT"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
cd "$ROOT"
while read -r p f; do
  [ -n "$p" ] || continue
  b="$(basename "$f" .R)"
  for arm in base lane; do
    "$RS" dev/superassign-trace.R "$arm" "$p" "$f" \
      "$OUT/$p--$b--$arm.tsv" > "$OUT/$p--$b--$arm.txt" 2>&1
  done
done < "$1"
echo "TRACE done"
