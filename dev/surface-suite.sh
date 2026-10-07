#!/usr/bin/env bash
# Lane surface: every test file of all eight packages, one R process
# per file, N at a time, on the lane build (dev/surface-runtest.R, every
# gate on). The extensions other than frmtmb.sample are the rellib-r6
# builds behind the lane library, so each loads the lane's frmtmb.
#
#   bash dev/surface-suite.sh <tag> [N] [pkg ...]
set -u
wt="$(cd "$(dirname "$0")/.." && pwd)"
tag="$1"; N="${2:-20}"; shift; shift || true
pkgs="${*:-frmtmb frmtmb.sample frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.spline}"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
: "${TMP:?TMP is unset}"
out="$wt/dev/surface-suite-$tag"
mkdir -p "$out"
for p in $pkgs; do
  if [ "$p" = frmtmb ]; then d="$wt/tests/testthat"; else d="$wt/extensions/$p/tests/testthat"; fi
  for f in "$d"/test-*.R; do
    while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
    "$RS" "$wt/dev/surface-runtest.R" "$p" "$f" lane > "$out/$p--$(basename "$f" .R).txt" 2>&1 &
  done
done
wait
echo "surface-suite $tag done: $(ls "$out" | wc -l) files"
