#!/usr/bin/env bash
# Lane silent2: one session of one suite, in one order, with the
# environment of dev/silent2-sessions.sh and run-sessions.sh, into
# dev/release/silent2-sess-sessions/<package>--<order>.txt.
#   bash dev/silent2-session1.sh <package> <order>
set -u
: "${TMP:?TMP is unset}"
p="$1"; o="$2"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-silent2-lib
unset R_PARALLEL_PORT
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
if [ "$p" = frmtmb ]; then d="$ROOT/tests/testthat"
else d="$ROOT/extensions/$p/tests/testthat"; fi
OUT="$ROOT/dev/release/silent2-sess-sessions"
mkdir -p "$OUT"
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" \
  "$ROOT/dev/release/run-session.R" "$p" "$o" "$d" > "$OUT/$p--$o.txt" 2>&1
