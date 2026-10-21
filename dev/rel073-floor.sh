#!/usr/bin/env bash
# Does frmtmb.sample 0.21.0 need frmtmb 0.73.0? It is installed ALONE
# into a scratch library over rellib-r10 (frmtmb 0.72.0 and the 0.72.0
# extensions), from a copy whose DESCRIPTION floor is lowered so that
# the install does not refuse on the version, and loaded there. If it
# installs and loads, every test file of its suite runs there, one per
# process, every gate on, as dev/rel072-floor.sh did; a file that fails
# only on the old core names the floor. Adapted from dev/rel072-floor.sh.
#
#   bash dev/rel073-floor.sh [N]
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/3}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/3}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
FLD=C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel073-log/floor-lib
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
cd "$ROOT"
p=frmtmb.sample
FL=$FLD/$p
rm -rf "$FL"; mkdir -p "$FL"
SRC=$FLD-src/$p; rm -rf "$SRC"; mkdir -p "$FLD-src"
cp -r "extensions/$p" "$SRC"
sed -i -e 's/frmtmb (>= 0.73.0)/frmtmb (>= 0.61.0)/' "$SRC/DESCRIPTION"
grep -n 'frmtmb (>=' "$SRC/DESCRIPTION"
echo "== install $p 0.21.0 over rellib-r10"
R_LIBS="$FL;C:/Users/adf44/source/r/rellib-r10;C:/Users/adf44/AppData/Local/R/win-library/4.6" \
  "$R" CMD INSTALL --library="$FL" --no-multiarch "$SRC" 2>&1 |
  grep -E 'DONE|ERROR|Error|error|object .* not exported|is not an exported|installation of package'
echo "== load"
"$RS" -e ".libPaths(c('$FL', 'C:/Users/adf44/source/r/rellib-r10', 'C:/Users/adf44/AppData/Local/R/win-library/4.6')); cat('frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ library(frmtmb.sample); paste('loaded', format(packageVersion('frmtmb.sample')), dirname(find.package('frmtmb.sample'))) }, error = function(e) paste('LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
if [ -d "$FL/$p" ]; then
  : > "dev/rel073-log/floor-$p.jobs"
  for f in extensions/$p/tests/testthat/test-*.R; do
    echo "$p $f" >> "dev/rel073-log/floor-$p.jobs"
  done
  ROUND_BASE=C:/Users/adf44/source/r/rellib-r10 /usr/bin/bash dev/ciharden-scan.sh \
    "r073-floor-$p" ref "$FL" "${1:-8}" "dev/rel073-log/floor-$p.jobs"
else
  echo "NOT INSTALLED: $p 0.21.0 does not install over frmtmb 0.72.0"
fi
