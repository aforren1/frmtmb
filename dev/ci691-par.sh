#!/usr/bin/env bash
# Lane ci691: run a list of test files, one R process per file, N at a
# time, through dev/ci691-run1.R against a chosen library, and write one
# summary line per file in list order.
#
#   bash dev/ci691-par.sh <tag> <job list> <lib|base> [N] [Rscript]
#
# The job list has one "<package> <test file>" line per file. Output:
# dev/ci691-suite-<tag>/<package>--<file>.txt per file, and
# dev/ci691-suite-<tag>.sum with the RESULT lines (gitignored).
set -u
tag="$1"; list="$2"; lib="$3"; N="${4:-16}"
RS="${5:-/c/Program Files/R/R-4.6.1/bin/Rscript.exe}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RUN="$ROOT/dev/ci691-run1.R"
OUT="$ROOT/dev/ci691-suite-$tag"
SUM="$ROOT/dev/ci691-suite-$tag/ci691-$tag.sum"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
rm -rf "$OUT"; mkdir -p "$OUT"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$RUN" "$lib" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
: > "$SUM"
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -n 1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$p $hit" >> "$SUM"
  else
    echo "$p NO RESULT LINE for $f" >> "$SUM"
  fi
done < "$list"
echo "SUITE $tag ran $ran of $want" >> "$SUM"
tail -n 1 "$SUM"
