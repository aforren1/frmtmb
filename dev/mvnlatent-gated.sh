#!/usr/bin/env bash
# The gated tier (FRMTMB_BRMS_FIT_TESTS=true) on the files that skip in
# the plain tier, one file per process, N at a time, through
# dev/mvnlatent-run1.R, on the lane's build.
#
#   bash dev/mvnlatent-gated.sh <list of "pkg file"> <tier> <N>
set -u
list="$1"; tier="$2"; N="$3"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/mvnlatent-out/$tier"; LOG="$OUT.log"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$LOG"
while read -r p f; do
  if [ "$p" = frmtmb ]; then d="$ROOT/tests/testthat"
  else d="$ROOT/extensions/$p/tests/testthat"; fi
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  ( cd "$d" && "$RS" "$ROOT/dev/mvnlatent-run1.R" "$p" "$d/$f" \
      > "$OUT/$p--${f%.R}.txt" 2>&1 ) &
done < "$list"
wait
while read -r p f; do
  h=$(grep -a "^RESULT" "$OUT/$p--${f%.R}.txt" | tail -1)
  echo "$p ${h:-NO RESULT LINE for $f}" >> "$LOG"
done < "$list"
grep -c RESULT "$LOG"
