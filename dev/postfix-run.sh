#!/usr/bin/env bash
# Lane postfix: run a list of test files, one R process per file, N at a
# time, through dev/release/run-tests.R, with the lane library in front
# of rellib-r8 (FRMTMB_RELLIB) unless the arm is "before".
#
#   bash dev/postfix-run.sh <arm> <tag> <job list> [N]
#
# The job list has one "<package> <test file>" line per file. Output:
# dev/postfix-tests/<tag>/<package>--<file>.txt per file and
# dev/postfix-tests/<tag>.log with the RESULT lines in list order. The
# caller sets NOT_CRAN and FRMTMB_BRMS_FIT_TESTS.
set -u
arm="$1"; tag="$2"; list="$3"; N="${4:-6}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# POSTFIX_RS runs the files under another R, the OpenBLAS copy of
# dev/ciharden-openblas.sh
RS="${POSTFIX_RS:-/c/Program Files/R/R-4.6.1/bin/Rscript.exe}"
RUN="$ROOT/dev/release/run-tests.R"
OUT="$ROOT/dev/postfix-tests/$tag"
LOG="$ROOT/dev/postfix-tests/$tag.log"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
if [ "$arm" = "before" ]; then
  unset FRMTMB_RELLIB
else
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-postfix-lib
fi
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$ROOT/$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( cd "$ROOT" && "$RS" "$RUN" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
want=0; ran=0
while read -r p f; do
  [ -n "$p" ] || continue
  want=$((want + 1))
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  lib=$(grep -a "^lib:" "$o" | head -1)
  if [ -n "$hit" ]; then
    ran=$((ran + 1)); echo "$p $hit | $lib" >> "$LOG"
  else
    echo "$p NO RESULT LINE for $f" >> "$LOG"; tail -5 "$o" >> "$LOG"
  fi
done < "$list"
echo "TAG $tag arm $arm ran $ran of $want" >> "$LOG"
