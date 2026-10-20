#!/usr/bin/env bash
# Do the merged extensions need frmtmb 0.72.0? Each of the four whose
# floor this release raises (frmtmb.sample, frmtmb.eam, frmtmb.latent,
# frmtmb.coupling) is installed ALONE into its own scratch library over
# rellib-r9 (frmtmb 0.71.0 and the 0.71.0 extensions), and every test
# file of its suite runs there, one per process, every gate on. A file
# that fails only on the old core names the floor. Adapted from
# dev/rel071-floor.sh.
#
#   bash dev/rel072-floor.sh [N]
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/3}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/3}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
FLD=C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel072-log/floor-lib
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
cd "$ROOT"
: > dev/rel072-log/floor.jobs
for p in frmtmb.sample frmtmb.eam frmtmb.latent frmtmb.coupling; do
  FL=$FLD/$p
  rm -rf "$FL"; mkdir -p "$FL"
  # a copy whose DESCRIPTION floor is lowered, so that the install does
  # not refuse frmtmb 0.71.0 and the tests decide
  SRC=$FLD-src/$p; rm -rf "$SRC"; mkdir -p "$FLD-src"
  cp -r "extensions/$p" "$SRC"
  sed -i -e 's/frmtmb (>= 0.72.0)/frmtmb (>= 0.61.0)/' \
    -e 's/frmtmb.sample (>= 0.20.0)/frmtmb.sample/' "$SRC/DESCRIPTION"
  R_LIBS="$FL;C:/Users/adf44/source/r/rellib-r9;C:/Users/adf44/AppData/Local/R/win-library/4.6" \
    "$R" CMD INSTALL --library="$FL" --no-multiarch "$SRC" 2>&1 |
    grep -E 'DONE|ERROR'
  : > "dev/rel072-log/floor-$p.jobs"
  for f in extensions/$p/tests/testthat/test-*.R; do
    echo "$p $f" >> "dev/rel072-log/floor-$p.jobs"
  done
done
for p in frmtmb.sample frmtmb.eam frmtmb.latent frmtmb.coupling; do
  ROUND_BASE=C:/Users/adf44/source/r/rellib-r9 /usr/bin/bash dev/ciharden-scan.sh \
    "r072-floor-$p" ref "$FLD/$p" "${1:-8}" "dev/rel072-log/floor-$p.jobs"
done
