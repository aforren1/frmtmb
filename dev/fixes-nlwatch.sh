#!/usr/bin/env bash
# Lane fixes, punch round (from the reviewer's dev/fixes-rev-nlwatch.sh): every test file that mentions a nonlinear
# formula, with the guard wrapped (dev/fixes-nlwatch-run.R).
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-fixes
LIB=C:/Users/adf44/source/r/wt-fixes-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
: "${TMP:?TMP is unset}"
OUT="$ROOT/dev/fixes-nlwatch-log"
mkdir -p "$OUT"
N=16
cd "$ROOT"
for f in $(grep -l "nl = TRUE\|nlf(\|nl=TRUE" tests/testthat/test-*.R extensions/*/tests/testthat/test-*.R); do
  case "$f" in
    extensions/*) p=$(echo "$f" | cut -d/ -f2); dir="$ROOT/extensions/$p" ;;
    *) p=frmtmb; dir="$ROOT" ;;
  esac
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  b=$(basename "$f" .R)
  ( cd "$dir" && "$RS" "$ROOT/dev/fixes-nlwatch-run.R" "$LIB" "$p" \
      "$ROOT/$f" "$OUT/$p--$b.calls" > "$OUT/$p--$b.txt" 2>&1 ) &
done
wait
echo done
