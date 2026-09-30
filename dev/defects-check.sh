#!/bin/sh
# Lane wt-defects: R CMD build and R CMD check --as-cran, once per
# changed package, inside dev/defects-check/.
#   sh dev/defects-check.sh core|sample
cd "$(dirname "$0")/.." || exit 1
W="$(pwd -W)"
case "$1" in
  core) src="$W"; name=frmtmb ;;
  sample) src="$W/extensions/frmtmb.sample"; name=frmtmb.sample ;;
esac
out=dev/defects-check/$1
rm -rf "$out"; mkdir -p "$out"; cd "$out" || exit 1
RC="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-defects-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_LIBS_USER="$R_LIBS"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
"$RC" CMD build "$src" > build.log 2>&1 || { echo BUILD FAILED; tail build.log; exit 1; }
tb=$(ls ${name}_*.tar.gz)
"$RC" CMD check --as-cran "$tb" > check.log 2>&1
grep -E "^Status|ERROR|WARNING|NOTE" check.log
