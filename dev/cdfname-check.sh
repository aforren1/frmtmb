#!/usr/bin/env bash
# R CMD check --as-cran, one package at a time, built with vignettes,
# inside dev/cdfname-check/ (as dev/release/run-check.ps1, on the lane's
# library).   bash dev/cdfname-check.sh core frmtmb.eam ...
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-cdfname
OUT=$WT/dev/cdfname-check
LOG=$WT/dev/cdfname-check.log
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export R_LIBS="C:/Users/adf44/source/r/wt-cdfname-lib;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
mkdir -p $OUT
cd $OUT || exit 1
for p in "$@"; do
  if [ "$p" = core ]; then d=$WT; n=frmtmb; else d=$WT/extensions/$p; n=$p; fi
  echo "===== CHECK $n =====" >> $LOG
  rm -f ${n}_*.tar.gz
  "$R" CMD build "$d" >> $LOG 2>&1
  tgz=$(ls -t ${n}_*.tar.gz 2>/dev/null | head -1)
  if [ -z "$tgz" ]; then echo "NO TARBALL BUILT" >> $LOG; continue; fi
  "$R" CMD check --as-cran "$tgz" >> $LOG 2>&1
  echo "CHECKED $n: $(grep -a '^Status' $n.Rcheck/00check.log)" >> $LOG
done
echo "CHECK END" >> $LOG
