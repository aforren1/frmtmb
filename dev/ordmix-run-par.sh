#!/usr/bin/env bash
# dev/release/run-par.sh for lane ordmix: a list of test files, one R
# process per file, N at a time, through dev/ordmix-runtest.R.
#
#   bash dev/ordmix-run-par.sh <arm> <tier> <job list> [N]
#
# The job list has one "<package> <test file>" line per file. The tier
# environment is set here: NOT_CRAN and every gate, so that nothing is
# skipped for want of a switch. Output: dev/ordmix-suite-<tier>/ with
# one log per file and <tier>.log with the RESULT lines in list order.
set -u
arm="$1"; tier="$2"; list="$3"; N="${4:-16}"
ROOT="C:/Users/adf44/source/r/frmtmb-wt-ordmix"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/ordmix-runtest.R"
OUT="$ROOT/dev/ordmix-suite-$tier"
LOG="$OUT/$tier.log"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true
: "${TMP:?TMP is unset, and R cannot create its temporary directory}"
# the native toolchain reads TMPDIR first, and Git Bash sets it to an
# MSYS path the compiler cannot use
export TMPDIR="$TMP"
rm -rf "$OUT"; mkdir -p "$OUT"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" "$RUN" "$p" "$f" "$arm" gated > "$o" 2>&1 ) &
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
    echo "$p NO RESULT LINE for $f" >> "$LOG"
  fi
done < "$list"
echo "RAN $ran of $want" >> "$LOG"
echo "RAN $ran of $want"
