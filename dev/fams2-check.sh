#!/bin/bash
# Build and R CMD check --as-cran one package inside dev/fams2-check/.
#   bash dev/fams2-check.sh <path to package dir>
set -e
root=/c/Users/adf44/source/r/frmtmb-wt-fams2
pkgdir=$1
out=$root/dev/fams2-check
mkdir -p "$out"
cd "$out"
export TMP="C:\Users\adf44\AppData\Local\Temp\1" TEMP="C:\Users\adf44\AppData\Local\Temp\1"
export RSTUDIO_PANDOC="C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="C:/Users/adf44/source/r/wt-fams2-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=
pkg=$(grep "^Package:" "$pkgdir/DESCRIPTION" | sed "s/Package: *//")
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD build "$pkgdir" > build-$pkg.log 2>&1
tb=$(ls -t ${pkg}_*.tar.gz | head -1)
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD check --as-cran "$tb" > check-$pkg.log 2>&1 || true
grep -h "Status:" check-$pkg.log
