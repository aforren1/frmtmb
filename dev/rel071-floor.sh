#!/usr/bin/env bash
# Does the merged frmtmb.eam need frmtmb 0.71.0? Install it alone into a
# scratch library over rellib-r8 (frmtmb 0.70.0) and run every test file
# of its suite there, one per process, every gate on. A file that fails
# only on the old core names the floor.
#
#   bash dev/rel071-floor.sh
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/3}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/3}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
FL=C:/Users/adf44/source/r/rel071-floor-lib
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$FL;C:/Users/adf44/source/r/rellib-r8;C:/Users/adf44/AppData/Local/R/win-library/4.6"
rm -rf "$FL"; mkdir -p "$FL"
"$R" CMD INSTALL --library="$FL" --no-multiarch "$ROOT/extensions/frmtmb.eam" 2>&1 |
  grep -E 'DONE|ERROR'
cd "$ROOT"
: > dev/rel071-log/floor.jobs
for f in extensions/frmtmb.eam/tests/testthat/test-*.R; do
  echo "frmtmb.eam $f" >> dev/rel071-log/floor.jobs
done
TMP="$TMP" ROUND_BASE=C:/Users/adf44/source/r/rellib-r8 /usr/bin/bash dev/ciharden-scan.sh \
  r071-floor-eam ref "$FL" 8 dev/rel071-log/floor.jobs
