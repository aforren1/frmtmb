#!/usr/bin/env bash
# Lane arstep: one tier of test files, one file per R process, N at a time,
# on the lane library ahead of rellib-r18. Args: tier joblist N
export TMP="C:/Users/adf44/AppData/Local/Temp" TEMP="C:/Users/adf44/AppData/Local/Temp"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-arstep-lib
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/stan-cache
cd /c/Users/adf44/source/r/frmtmb-wt-arstep
bash dev/arstep-run-par.sh "$1" "$2" "$3"
