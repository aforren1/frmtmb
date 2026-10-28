#!/bin/bash
# One-session runs (forward, reverse, preload) on the lane library.
# Usage: bash dev/reslot-sessions.sh <tier> <N> [package ...]
export NOT_CRAN=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-reslot-lib
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/stan-cache
bash C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/release/run-sessions.sh "$@"
