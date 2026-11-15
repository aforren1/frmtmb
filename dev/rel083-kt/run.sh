#!/bin/bash
# usage: run.sh <workdir> <outfile> <script> [args...]; runs one Rscript in the foreground
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win NOT_CRAN=true
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/stan-cache
export MC_CORES=1 OMP_NUM_THREADS=1
wd=$1; out=$2; shift 2
cd "$wd" || exit 1
start=$(date +%s)
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$@" > "$out" 2>&1
rc=$?
echo "## exit $rc, $(( $(date +%s) - start )) s" >> "$out"
exit $rc
