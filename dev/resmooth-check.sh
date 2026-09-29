#!/bin/sh
# Lane wt-resmooth: the single R CMD check --as-cran of this round.
# pandoc and TinyTeX on PATH, or check reports a bogus pdflatex ERROR
# (dev/lane-rules.md). No --no-manual: the HTML and PDF manual sections
# are where an Rd defect surfaces.
cd /c/Users/adf44/source/r
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-resmooth-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
rm -rf /c/Users/adf44/source/r/resmooth-check
mkdir -p /c/Users/adf44/source/r/resmooth-check
cd /c/Users/adf44/source/r/resmooth-check
"C:/Program Files/R/R-4.6.1/bin/R.exe" CMD build \
  /c/Users/adf44/source/r/frmtmb-wt-resmooth
tb=$(ls frmtmb_*.tar.gz | head -1)
echo "BUILT $tb"
"C:/Program Files/R/R-4.6.1/bin/R.exe" CMD check --as-cran "$tb"
echo CHECK-DONE
