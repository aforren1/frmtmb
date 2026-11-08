#!/usr/bin/env bash
# R CMD build (with vignettes) and R CMD check --as-cran, once, for each
# package this lane changed, one after the other, inside
# dev/superassign-check/. Flags as dev/release/run-check.ps1.
#
#   bash dev/superassign-check.sh
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
OUT="$ROOT/dev/superassign-check"
export R_LIBS="C:/Users/adf44/source/r/wt-superassign-lib;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
for p in "$ROOT" "$ROOT/extensions/frmtmb.ode" "$ROOT/extensions/frmtmb.sample" \
         "$ROOT/extensions/frmtmb.survival"; do
  n=$(basename "$p"); [ "$p" = "$ROOT" ] && n=frmtmb
  echo "===== CHECK $n ====="
  "$R" CMD build "$p" > "build-$n.log" 2>&1
  tgz=$(ls -t ${n}_*.tar.gz 2>/dev/null | head -1)
  if [ -z "$tgz" ]; then echo "NO TARBALL $n"; tail -5 "build-$n.log"; continue; fi
  "$R" CMD check --as-cran "$tgz" > "check-$n.log" 2>&1
  grep -E '^Status:' "check-$n.log" | sed "s/^/$n /"
done
echo "CHECKS DONE"
