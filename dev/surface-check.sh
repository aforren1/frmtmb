#!/usr/bin/env bash
# Lane surface: R CMD build and check --as-cran of core and
# frmtmb.sample, one at a time, inside dev/surface-check/.
#
#   bash dev/surface-check.sh
set -u
wt="$(cd "$(dirname "$0")/.." && pwd -W)"
out="$wt/dev/surface-check"
mkdir -p "$out"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_LIBS="C:/Users/adf44/source/r/wt-surface-lib;C:/Users/adf44/source/r/rellib-r6;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export NOT_CRAN=false
cd "$out"
for p in ${CHECK_PKGS:-frmtmb frmtmb.sample}; do
  if [ "$p" = frmtmb ]; then src="$wt"; else src="$wt/extensions/$p"; fi
  rm -f ${p}_*.tar.gz
  # with vignettes, as dev/release/run-check.ps1 builds
  "$R" CMD build "$src" > "build-$p.log" 2>&1
  tb="$(ls ${p}_*.tar.gz | head -1)"
  "$R" CMD check --as-cran "$tb" > "check-$p.log" 2>&1
  grep -h "^Status" "$p.Rcheck/00check.log" || tail -5 "check-$p.log"
done
echo "surface-check done"
