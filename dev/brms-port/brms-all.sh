#!/usr/bin/env bash
# Fit every audited vignette's models with brms, one process per
# vignette, all at once. Compiles dominate, and they run in parallel.
#
#   PORT_OUT=<dir> dev/brms-port/brms-all.sh [iter] [chains]
set -u
here="$(cd "$(dirname "$0")" && pwd)"
iter="${1:-1000}"
chains="${2:-2}"
rscript="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
mkdir -p "$PORT_OUT/brms"
for v in brms_overview brms_multilevel brms_distreg brms_nonlinear \
         brms_phylogenetics brms_monotonic brms_multivariate \
         brms_missings brms_customfamilies; do
  "$rscript" "$here/brms-fit.R" "$v" "$iter" "$chains" \
    > "$PORT_OUT/brms/$v.out" 2>&1 &
done
wait
grep -h "^done:" "$PORT_OUT"/brms/*.log
