#!/usr/bin/env bash
# Install core, or a package under extensions/, into the lane's library.
#   bash dev/refuse-install.sh [. | extensions/frmtmb.sample]
src="${1:-.}"
export R_LIBS="C:/Users/adf44/source/r/wt-refuse-lib;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
( PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
  "/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --no-multiarch --no-test-load \
    --library="C:/Users/adf44/source/r/wt-refuse-lib" "$src" 2>&1 | tail -2 )
