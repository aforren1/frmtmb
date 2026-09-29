#!/bin/sh
# Reviewer: the gated brms likelihood file, fresh Stan cache, so every
# program in it is compiled here rather than served from the worker's.
export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$1"
mkdir -p "$FRMTMB_STAN_CACHE"
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" \
  C:/Users/adf44/source/r/frmtmb-wt-csfactor/dev/csfactor-rev-run.R \
  C:/Users/adf44/source/r/wt-csfactor-lib test-brms-likelihood.R
