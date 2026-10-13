#!/usr/bin/env bash
# Run dev/moparity-audit-run.R over every test file of the eight
# packages, one R process per file, N at a time.
#   bash dev/moparity-audit.sh [N]
N="${1:-16}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-moparity
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$ROOT/dev/moparity-audit-out
rm -rf "$OUT"; mkdir -p "$OUT"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
for pkg in frmtmb frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline; do
  if [ "$pkg" = frmtmb ]; then d=$ROOT/tests/testthat; else d=$ROOT/extensions/$pkg/tests/testthat; fi
  for f in "$d"/test-*.R; do
    while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
    b=$(basename "$f" .R)
    ( cd "$d" && "$RS" "$ROOT/dev/moparity-audit-run.R" "$pkg" "$f" "$OUT/$pkg--$b.res" > "$OUT/$pkg--$b.log" 2>&1 ) &
  done
done
wait
cat "$OUT"/*.res > "$OUT/all.txt"
echo "files with a result: $(wc -l < "$OUT/all.txt")"
