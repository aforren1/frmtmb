#!/usr/bin/env bash
# Lane surface: the mechanical port (raw and spell passes) and the hand
# translation against one build.
#
#   bash dev/surface-port.sh <tag> <libs>
#
# <tag> names the output directories dev/surface-port-out/<tag> and
# dev/surface-bv-out/<tag>; <libs> is the ";"-separated library list
# put before the user library (PORT_LIB and BV_LIB).
set -u
wt="$(cd "$(dirname "$0")/.." && pwd)"
tag="$1"
lib="$2"
rs="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PORT_LIB="$lib"
export PORT_OUT="$wt/dev/surface-port-out/$tag"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="$wt/dev/stan-cache"
export NOT_CRAN=true
mkdir -p "$PORT_OUT"
"$rs" "$wt/dev/brms-port/run-all.R" 120 raw > "$PORT_OUT/run-raw.log" 2>&1
"$rs" "$wt/dev/brms-port/run-all.R" 120 spell > "$PORT_OUT/run-spell.log" 2>&1
"$rs" "$wt/dev/brms-port/summarize.R" > "$PORT_OUT/summarize.log" 2>&1
bash "$wt/dev/brms-vignettes/_run-all.sh" "$lib" \
  "$wt/dev/surface-bv-out/$tag" > "$wt/dev/surface-port-out/$tag-bv.log" 2>&1
echo "surface-port $tag done"
