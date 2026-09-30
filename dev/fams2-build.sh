#!/bin/sh
# Roxygenise and install core (and frmtmb.sample when asked) into the lane library.
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
cd /c/Users/adf44/source/r/frmtmb-wt-fams2
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" -e '.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib","C:/Users/adf44/source/r/rellib-r3","C:/Users/adf44/AppData/Local/R/win-library/4.6")); roxygen2::roxygenise(".")' > dev/fams2-roxy.log 2>&1 || { tail -20 dev/fams2-roxy.log; exit 1; }
grep -i "warn\|error" dev/fams2-roxy.log | head
cd ..
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/wt-fams2-lib --no-multiarch frmtmb-wt-fams2 > frmtmb-wt-fams2/dev/fams2-install.log 2>&1
tail -1 frmtmb-wt-fams2/dev/fams2-install.log
if [ "$1" = "sample" ]; then
  cd frmtmb-wt-fams2/extensions/frmtmb.sample
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" -e '.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib","C:/Users/adf44/source/r/rellib-r3","C:/Users/adf44/AppData/Local/R/win-library/4.6")); roxygen2::roxygenise(".")' > ../../dev/fams2-roxy-sample.log 2>&1
  cd ..
  "/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/wt-fams2-lib --no-multiarch frmtmb.sample > ../dev/fams2-install-sample.log 2>&1
  tail -1 ../dev/fams2-install-sample.log
fi
