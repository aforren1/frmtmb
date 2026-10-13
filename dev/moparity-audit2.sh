#!/usr/bin/env bash
# dev/moparity-audit-run2.R over the files whose smallest theta of the
# first audit fell below -308.45 (the list on stdin as "pkg file").
#   bash dev/moparity-audit2.sh < list
ROOT=/c/Users/adf44/source/r/frmtmb-wt-moparity
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$ROOT/dev/moparity-audit2-out
rm -rf "$OUT"; mkdir -p "$OUT"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
while read -r pkg f; do
  [ -n "$pkg" ] || continue
  if [ "$pkg" = frmtmb ]; then d=$ROOT/tests/testthat; else d=$ROOT/extensions/$pkg/tests/testthat; fi
  b=$(basename "$f" .R)
  ( cd "$d" && "$RS" "$ROOT/dev/moparity-audit-run2.R" "$pkg" "$d/$f" "$OUT/$pkg--$b.res" > "$OUT/$pkg--$b.log" 2>&1 ) &
done
wait
cat "$OUT"/*.res > "$OUT/all.txt"
cat "$OUT/all.txt"
