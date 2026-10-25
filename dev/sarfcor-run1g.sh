#!/bin/sh
# Gated run of one test file: sh dev/sarfcor-run1g.sh <pkg> <file> <log>
export FRMTMB_BRMS_FIT_TESTS=true
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-sarfcor/dev/stan-cache
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" C:/Users/adf44/source/r/frmtmb-wt-sarfcor/dev/sarfcor-run1.R "$1" "$2" > "$3" 2>&1
