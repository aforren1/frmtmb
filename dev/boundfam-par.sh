#!/usr/bin/env bash
# Lane boundfam: run a job list through dev/boundfam-run-traced.R, three
# R processes at a time, one test file per process.
#
#   bash dev/boundfam-par.sh <arm> <job list> [N]
#
# <arm> "base" runs rellib-r20's frmtmb, "lane" the lane library's. Each
# file writes dev/boundfam-suite-<arm>/<pkg>--<file>.txt (the run-tests.R
# log) and .trace (dev/boundfam-run-traced.R's record). The caller sets
# NOT_CRAN and any gate.
set -u
arm="$1"; list="$2"; N="${3:-3}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="$ROOT/dev/boundfam-run-traced.R"
OUT="$ROOT/dev/boundfam-suite-$arm"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
: "${TMP:?TMP is unset}"
if [ "$arm" = "lane" ]; then
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-boundfam-lib
else
  export FRMTMB_RELLIB=
fi
mkdir -p "$OUT"
port=21800
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; exit 1; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  b="$OUT/$p--$(basename "$f" .R)"
  port=$((port + 1))
  ( BOUNDFAM_TRACE="$b.trace" R_PARALLEL_PORT=$port \
      "$RS" "$RUN" "$p" "$f" > "$b.txt" 2>&1 ) &
done < "$list"
wait
echo "done $arm"
