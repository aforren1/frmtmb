#!/usr/bin/env bash
# Reviewer of lane fixes: run named test files on a library, one R
# process per file, gated, with a COPY of the stan cache.
#
#   bash dev/fixes-rev-tests.sh <lib> <tag> <pkg> <file> [<pkg> <file>]...
set -u
LIB="$1"; TAG="$2"; shift 2
ROOT=/c/Users/adf44/source/r/frmtmb-wt-fixes
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/fixes-rev-run-test.R"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/fixes-rev-stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
: "${TMP:?TMP is unset}"
OUT="$ROOT/dev/fixes-rev-log"
mkdir -p "$OUT"
while [ $# -ge 2 ]; do
  p="$1"; f="$2"; shift 2
  if [ "$p" = frmtmb ]; then dir="$ROOT"; else dir="$ROOT/extensions/$p"; fi
  o="$OUT/$TAG--$p--$(basename "$f" .R).txt"
  ( cd "$dir" && "$RS" "$RUN" "$LIB" "$p" "$dir/tests/testthat/$f" > "$o" 2>&1 ) &
done
wait
grep -a -H "^RESULT\|^lib:" "$OUT/$TAG--"*.txt
