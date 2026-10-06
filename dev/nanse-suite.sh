#!/usr/bin/env bash
# All eight suites, one R process per test file, N at a time, against
# the lane library, every gate on. Per file: <out>/<pkg>--<file>.txt
# (the run) and <out>/<pkg>--<file>.fire (one line per firing of the
# standard-error warning). Summary: <out>/suite.log.
#   bash dev/nanse-suite.sh <out dir> [N] [pkg ...]
set -u
OUT="$1"; N="${2:-20}"; shift 2 || true
W=/c/Users/adf44/source/r/frmtmb-wt-nanse
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$W/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true
: "${TMP:?TMP is unset}"
pkgs="${*:-frmtmb frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline}"
mkdir -p "$W/$OUT"
LIST="$W/$OUT/jobs.txt"; : > "$LIST"
for p in $pkgs; do
  if [ "$p" = frmtmb ]; then d="$W/tests/testthat"; else
    d="$W/extensions/$p/tests/testthat"; fi
  for f in "$d"/test-*.R; do echo "$p $f" >> "$LIST"; done
done
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  b="$p--$(basename "$f" .R)"
  ( cd "$(dirname "$f")" && "$RS" "$W/dev/nanse-run-tests.R" "$p" "$f" \
      "$W/$OUT/$b.fire" > "$W/$OUT/$b.txt" 2>&1 ) &
done < "$LIST"
wait
LOG="$W/$OUT/suite.log"; : > "$LOG"
while read -r p f; do
  b="$p--$(basename "$f" .R)"
  hit=$(grep -a "^RESULT" "$W/$OUT/$b.txt" | tail -1)
  if [ -n "$hit" ]; then echo "$p $hit" >> "$LOG"; else
    echo "$p NO RESULT LINE for $f" >> "$LOG"; fi
done < "$LIST"
echo "ran $(grep -c RESULT "$LOG") of $(wc -l < "$LIST")" >> "$LOG"
