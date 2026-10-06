#!/usr/bin/env bash
# Reviewer: rerun the mechanical port on rellib-r5 with the lane's
# unchanged runner (seeded) and with the reviewer's unseeded copy.
#
#   bash dev/vigport-rev-mech.sh
set -u
dev="C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
rs="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PORT_LIB="C:/Users/adf44/source/r/rellib-r5"
export FRMTMB_STAN_CACHE="$dev/vigport-rev-out/stan-cache"
export PORT_PATCHES="brms_nonlinear.9.1,brms_multilevel.15.3"
VIGS="brms_overview brms_multilevel brms_distreg brms_nonlinear
  brms_phylogenetics brms_monotonic brms_multivariate brms_missings
  brms_customfamilies"
run_arm() {
  arm="$1"; script="$2"; shift 2
  export PORT_OUT="$dev/vigport-rev-out/$arm"
  mkdir -p "$PORT_OUT"
  for mode in "$@"; do
    for v in $VIGS; do
      "$rs" "$script" "$v" 120 "$mode" > "$PORT_OUT/$v-$mode.out" 2>&1 &
    done
    wait
  done
  "$rs" "$dev/brms-port/summarize.R" > "$PORT_OUT/summary.txt" 2>&1
}
run_arm seeded "$dev/brms-port/run-vignette.R" raw spell need keepprior &
run_arm noseed "$dev/vigport-rev-out/run-vignette-noseed.R" raw spell &
wait
echo done
