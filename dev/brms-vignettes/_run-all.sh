#!/usr/bin/env bash
# Run every hand-translated vignette script against one installed build,
# one R process per script, all at once.
#
#   dev/brms-vignettes/_run-all.sh <libs> <outdir>
#
# <libs> is BV_LIB (";"-separated libraries before the user library);
# <outdir> receives one CSV per script (BV_OUT) and one log per script.
# Then run _scoreboard.R and _summarize.R with BV_OUT=<outdir>.
#
# Why all at once: the scripts share nothing but the read-only library,
# and the slowest one decides the wall time either way.
set -u
here="$(cd "$(dirname "$0")" && pwd -W 2>/dev/null || pwd)"
lib="$1"
out="$2"
mkdir -p "$out"
rscript="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export BV_LIB="$lib" BV_OUT="$out" BV_DIR="$here"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="${FRMTMB_STAN_CACHE:-$here/../stan-cache}"
export NOT_CRAN=true
for f in "$here"/brms_*.R "$here"/brms-*.R; do
  nm="$(basename "$f" .R)"
  "$rscript" "$f" > "$out/bv-$nm.log" 2>&1 &
done
wait
for f in "$out"/bv-*.log; do
  ok="$(grep -c '^\[ok' "$f")"
  err="$(grep -c '^\[err' "$f")"
  done_line="$(grep -m1 '^### done' "$f" || echo 'NO DONE LINE')"
  printf '%-40s %s ok, %s err, %s\n' "$(basename "$f")" "$ok" "$err" \
    "$done_line"
done
