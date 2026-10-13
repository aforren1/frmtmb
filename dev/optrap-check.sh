#!/usr/bin/env bash
# Lane optrap: R CMD build and check --as-cran of core, inside
# dev/optrap-check/, against the lane's library.
#   bash dev/optrap-check.sh
set -eu
WT=C:/Users/adf44/source/r/frmtmb-wt-optrap
OUT=/c/Users/adf44/source/r/frmtmb-wt-optrap/dev/optrap-check
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="C:/Users/adf44/source/r/wt-optrap-lib;C:/Users/adf44/source/r/rellib-r7;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
# a stray test artifact in the tree would be built into the tarball
for f in tests/testthat/testthat-problems.rds tests/testthat/_problems \
         Rplots.pdf tests/testthat/Rplots.pdf; do
  if [ -e "$WT/$f" ]; then echo "stray $f; remove it first"; exit 1; fi
done
mkdir -p "$OUT"
cd "$OUT"
"$R" CMD build "$WT" > build.log 2>&1
tgz=$(ls -t frmtmb_*.tar.gz | head -1)
"$R" CMD check --as-cran "$tgz" > check.log 2>&1 || true
grep -E "^Status|^\* checking .* \.\.\. (WARNING|ERROR|NOTE)" check.log || true
tail -3 check.log
