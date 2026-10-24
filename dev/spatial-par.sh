#!/usr/bin/env bash
# Run a job list ("<package> <test file>" per line), one R process per
# file, N at a time, through dev/spatial-run-tests.R.
#
#   SPATIAL_LIBS="libA;libB" bash dev/spatial-par.sh <tag> <job list> [N]
#
# Writes dev/spatial-log/<tag>/<package>--<file>.txt and
# dev/spatial-log/<tag>.log (the RESULT lines in list order). The per-file
# logs are written by each R process's own redirection; nothing tails
# them.
set -u
tag="$1"; list="$2"; N="${3:-8}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/spatial-run-tests.R"
OUT="$ROOT/dev/spatial-log/$tag"
LOG="$ROOT/dev/spatial-log/$tag.log"
: "${SPATIAL_LIBS:?set SPATIAL_LIBS}"
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER="${R_MAKEVARS_USER:-C:/Users/adf44/Documents/.R/Makevars.win}"
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$ROOT/dev/stan-cache}"
export SPATIAL_LIBS
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$RUN" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$p $hit" >> "$LOG"
  else
    echo "$p NO RESULT LINE for $f" >> "$LOG"; tail -5 "$o" >> "$LOG"
  fi
done < "$list"
echo "RAN $ran of $want" >> "$LOG"
echo "RAN $ran of $want"
