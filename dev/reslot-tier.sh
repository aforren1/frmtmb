#!/bin/bash
# Run a job list on the lane library (FRMTMB_RELLIB) or on rellib-r14.
# Usage: bash dev/reslot-tier.sh <tier> <jobs> <after|before> [N] [gated]
tier=$1; jobs=$2; arm=$3; N=${4:-6}
export NOT_CRAN=true
if [ "${5:-}" = gated ]; then export FRMTMB_BRMS_FIT_TESTS=true; fi
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/stan-cache
if [ "$arm" = after ]; then
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-reslot-lib
else
  unset FRMTMB_RELLIB
fi
bash C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/release/run-par.sh "$tier" "$jobs" "$N"
