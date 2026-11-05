#!/usr/bin/env bash
# Item 3: an independent rerun, absolute paths, the five changed pairs and ten
# others (seed in the choice: every 16th of the sample), into my own directory.
export FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-refuse-lib;C:/Users/adf44/source/r/rellib-r16"
export FRMTMB_STAN_CACHE="$(cd dev/review-refuse-out/stan-cache && pwd -W)"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win" NOT_CRAN=true
O=dev/review-refuse-f2-equiv-out; mkdir -p $O/logs
for id in p029 p043 p070 p087 p169 p005 p021 p037 p053 p069 p085 p101 p117 p133 p149; do echo $id; done | xargs -P 2 -I{} bash -c '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/equiv-measure.R {} dev/review-refuse-f2-equiv-out --force > dev/review-refuse-f2-equiv-out/logs/{}.log 2>&1'
touch $O/DONE
