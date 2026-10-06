#!/usr/bin/env bash
# Reviewer, re-check: rerun the mechanical spell pass on rellib-r5 with
# the lane's punch-round patches.R, into the reviewer's own directory.
#
#   bash dev/vigport-rev2-spell.sh
set -u
dev="C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
rs="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PORT_LIB="C:/Users/adf44/source/r/rellib-r5"
export FRMTMB_STAN_CACHE="$dev/vigport-rev2-out/stan-cache"
export PORT_OUT="$dev/vigport-rev2-out/spell"
mkdir -p "$PORT_OUT" "$FRMTMB_STAN_CACHE"
for v in brms_overview brms_multilevel brms_distreg brms_nonlinear \
         brms_phylogenetics brms_monotonic brms_multivariate \
         brms_missings brms_customfamilies; do
  "$rs" "$dev/brms-port/run-vignette.R" "$v" 120 spell \
    > "$PORT_OUT/$v.out" 2>&1 &
done
wait
echo done
