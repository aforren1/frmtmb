#!/usr/bin/env bash
# Reviewer of lane fixes (re-check): every test file of core and the seven extensions, one R
# process per file, N at a time, against a chosen library.
#
#   bash dev/fixes-run-all.sh <lib> <outdir> [N] [list]
#
# Without a list it runs every tests/testthat/test-*.R of core and of
# each extension. Each process writes its own log; the summary is the
# RESULT line of each, in list order, in <outdir>/summary.txt.
set -u
LIB="$1"; OUT="$2"; N="${3:-16}"; LIST="${4:-}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/fixes-rev-run-test.R"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/fixes-rev2-stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
if [ -z "$LIST" ]; then
  LIST="$OUT/jobs.txt"
  : > "$LIST"
  for f in "$ROOT"/tests/testthat/test-*.R; do
    echo "frmtmb $ROOT $f" >> "$LIST"
  done
  for d in "$ROOT"/extensions/*/; do
    p="$(basename "$d")"
    for f in "$d"tests/testthat/test-*.R; do
      echo "$p $d $f" >> "$LIST"
    done
  done
fi
while read -r p dir f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( cd "$dir" && "$RS" "$RUN" "$LIB" "$p" "$f" > "$o" 2>&1 ) &
done < "$LIST"
wait
: > "$OUT/summary.txt"
while read -r p dir f; do
  [ -n "$p" ] || continue
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then
    echo "$p $hit" >> "$OUT/summary.txt"
  else
    echo "$p NO RESULT LINE for $(basename "$f")" >> "$OUT/summary.txt"
  fi
done < "$LIST"
echo "ran $(grep -c RESULT "$OUT/summary.txt") of $(grep -c . "$LIST")"
