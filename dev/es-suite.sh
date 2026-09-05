#!/bin/sh
# full core suite, one file per process, counts audited by name
ES_LIB="C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/es-lib"
WT="C:/Users/adf44/source/r/frmtmb-wt-esicar"
for f in $(ls "$WT"/tests/testthat/test-*.R | xargs -n1 basename); do
  R_LIBS="$ES_LIB;C:/Users/adf44/AppData/Local/R/win-library/4.6" NOT_CRAN=true \
    "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" "$WT/dev/es-runtest.R" "$f" 2>&1 |
    grep -E "^FILE|BAD"
done
