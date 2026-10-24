#!/usr/bin/env bash
# One R process per test file, N at a time, through
# dev/mvnlatent-run1.R, which also prints every test_that() block's
# time. A copy of dev/release/run-par.sh's loop with a test root and a
# package list as arguments, so the "before" arm can read a snapshot of
# the base tree.
#
#   bash dev/mvnlatent-par.sh <tier> <test root> <N> <package> ...
#
# The caller sets NOT_CRAN and the gates, and FRMTMB_RELLIB for the
# lane's build. Output: dev/mvnlatent-out/<tier>/<pkg>--<file>.txt and
# dev/mvnlatent-out/<tier>.log with one RESULT line per file.
set -u
tier="$1"; TR="$2"; N="$3"; shift 3
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/mvnlatent-run1.R"
OUT="$ROOT/dev/mvnlatent-out/$tier"
LOG="$ROOT/dev/mvnlatent-out/$tier.log"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER="${R_MAKEVARS_USER:-C:/Users/adf44/Documents/.R/Makevars.win}"
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$ROOT/dev/stan-cache}"
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
list="$OUT/jobs.txt"; : > "$list"
for p in "$@"; do
  if [ "$p" = frmtmb ]; then d="$TR/tests/testthat"
  else d="$TR/extensions/$p/tests/testthat"; fi
  for f in "$d"/test-*.R; do echo "$p $f" >> "$list"; done
done
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( cd "$(dirname "$f")" && "$RS" "$RUN" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r p f; do
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$p $hit" >> "$LOG"
  else
    echo "$p NO RESULT LINE for $f" >> "$LOG"
  fi
done < "$list"
echo "ran $ran of $want" >> "$LOG"
echo "ran $ran of $want"
