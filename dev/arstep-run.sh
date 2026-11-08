#!/usr/bin/env bash
# Lane arstep: run one R script with the round's environment. Args: log script args...
log=$1; shift
export TMP="C:/Users/adf44/AppData/Local/Temp" TEMP="C:/Users/adf44/AppData/Local/Temp"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/stan-cache
cd /c/Users/adf44/source/r/frmtmb-wt-arstep
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$@" > "$log" 2>&1
