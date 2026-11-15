#!/bin/bash
# One traced test file at a time (sequential), merged build.
R=/c/Users/adf44/source/r/frmtmb-wt-rel083
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win NOT_CRAN=true
export FRMTMB_STAN_CACHE=$R/dev/stan-cache MC_CORES=1 OMP_NUM_THREADS=1
unset FRMTMB_RELLIB FRMTMB_BRMS_FIT_TESTS
cd $R
OUT=C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-log/kt/selimit2/suite-r22
while read -r p f; do
  b="$OUT/$p--$(basename "$f" .R)"
  SELIMIT2_TRACE="$b.trace" R_PARALLEL_PORT=23901 "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/rel083-kt/selimit2/selimit2-run-traced.R "$p" "$f" > "$b.txt" 2>&1
  echo "$(date +%T) $p $f exit $?"
done < dev/rel083-kt/selimit2/jobs33.txt
echo ALLDONE
