#!/usr/bin/env bash
# Lane setier: test files one R process each, N at a time, every gate on
# (the scale measurements excepted, as in the release's tiers).
#   bash dev/setier-suite.sh <run name> <ref|ob> <libs> <job list|all> [N]
# <libs> is ';'-separated; dev/setier-run1.R appends the user library.
# A job list has "<package> <test file>" lines; "all" builds the eight
# suites from this worktree. Per file: dev/setier-log/<run>/<pkg>--<file>.txt
# and .cond (every frmtmb condition raised). Summary: <run>.sum beside it.
set -u
run="$1"; kind="$2"; libs="$3"; list="$4"; N="${5:-12}"
W=/c/Users/adf44/source/r/frmtmb-wt-setier
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$W/dev/stan-cache"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true
: "${TMP:?TMP is unset}"
if [ "$kind" = ob ]; then
  RS="$W/dev/setier-out/Rob/bin/x64/Rscript.exe"
  export OPENBLAS_NUM_THREADS=4
else
  RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
fi
OUT="$W/dev/setier-log/$run"
SUM="$W/dev/setier-log/$run.sum"
rm -rf "$OUT"; mkdir -p "$OUT"; rm -f "$SUM"
if [ "$list" = all ]; then
  list="$OUT/jobs.txt"; : > "$list"
  for p in frmtmb frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn \
           frmtmb.ode frmtmb.sample frmtmb.spline; do
    if [ "$p" = frmtmb ]; then d="$W/tests/testthat"; else
      d="$W/extensions/$p/tests/testthat"; fi
    for f in "$d"/test-*.R; do echo "$p $f" >> "$list"; done
  done
fi
RUN=$(cygpath -w "$W/dev/setier-run1.R")
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  b="$p--$(basename "$f" .R)"
  ( cd "$(dirname "$f")" && "$RS" "$RUN" "$libs" "$p" "$(cygpath -w "$f")" \
      "$(cygpath -w "$OUT/$b.cond")" > "$OUT/$b.txt" 2>&1 ) &
done < "$list"
wait
while read -r p f; do
  [ -n "$p" ] || continue
  o="$OUT/$p--$(basename "$f" .R).txt"
  hit=$(grep -a "^RESULT" "$o" | tail -1)
  if [ -n "$hit" ]; then echo "$p $hit" >> "$SUM"; else
    echo "$p NO RESULT LINE for $f" >> "$SUM"; tail -3 "$o" >> "$SUM"; fi
done < "$list"
echo "ran $(grep -c ' RESULT ' "$SUM") of $(grep -c . "$list")" >> "$SUM"
