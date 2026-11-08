#!/usr/bin/env bash
# Lane osa: roxygenise core and install it (and, with "ext", every
# extension) into the lane's private library.
#   bash dev/osa-install.sh [ext]
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
LIB="C:/Users/adf44/source/r/wt-osa-lib"
WT="C:/Users/adf44/source/r/frmtmb-wt-osa"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
R="/c/Program Files/R/R-4.6.1/bin"
"$R/Rscript.exe" -e "roxygen2::roxygenise('$WT')" || exit 1
"$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT" || exit 1
if [ "$1" = "ext" ]; then
  for p in coupling eam latent learn ode sample spatial survival; do
    "$R/Rscript.exe" -e "roxygen2::roxygenise('$WT/extensions/frmtmb.$p')" || exit 1
    "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/frmtmb.$p" || exit 1
  done
fi
