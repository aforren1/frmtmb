#!/usr/bin/env bash
# Lane nlband: run a tier's job list through dev/release/run-par.sh on the
# lane's library. Usage: bash dev/nlband-tier.sh <tier> <jobs> <N> [gated]
cd /c/Users/adf44/source/r/frmtmb-wt-nlband || exit 1
export NOT_CRAN=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-nlband-lib
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-nlband/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
if [ "${4:-}" = gated ]; then
  export FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true
fi
bash dev/release/run-par.sh "$1" "$2" "$3"
