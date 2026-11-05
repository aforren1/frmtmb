#!/usr/bin/env bash
# Short check of fix round 2, item 9: R CMD check --as-cran on core, frmtmb.sample and
# frmtmb.latent, one at a time, inside dev/review-predread-check/.
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-predread
mkdir -p "$W/dev/review-predread-check"
cd "$W/dev/review-predread-check" || exit 1
RB="/c/Program Files/R/R-4.6.1/bin/R.exe"
PD="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PD:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-predread-lib;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
O="$W/dev/review-predread-out/r2/rcheck"
mkdir -p "$O"
one() {
  nm="$1"; path="$2"
  "$RB" CMD build "$path" > "$O/$nm-build.txt" 2>&1
  tb=$(ls -t "$nm"_*.tar.gz | head -1)
  "$RB" CMD check --as-cran --no-multiarch "$tb" > "$O/$nm-check.txt" 2>&1
  echo "$nm finished" >> "$O/done.txt"
}
one frmtmb "$W"
one frmtmb.sample "$W/extensions/frmtmb.sample"
one frmtmb.latent "$W/extensions/frmtmb.latent"
echo finished > "$O/all.done"
