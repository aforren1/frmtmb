#!/usr/bin/env bash
# Review of lane boundfam: the gated brms tier, one test file per R
# process, through dev/release/run-tests.R.
#   bash dev/reviews/boundfam-review-gated.sh <lane|base> <job list> <outdir> [N] [port]
set -u
arm="$1"; list="$2"; OUT="$3"; N="${4:-2}"; port="${5:-23500}"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RUN="C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/release/run-tests.R"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
if [ "$arm" = "lane" ]; then
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-boundfam-lib
else
  export FRMTMB_RELLIB=
fi
mkdir -p "$OUT"
while read -r p f; do
  [ -n "$p" ] || continue
  [ -f "$f" ] || { echo "test file missing: $f"; continue; }
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  b="$OUT/$p--$(basename "$f" .R)"
  port=$((port + 1))
  ( R_PARALLEL_PORT=$port "$RS" "$RUN" "$p" "$f" > "$b.txt" 2>&1 ) &
done < "$list"
wait
echo "done $arm"
