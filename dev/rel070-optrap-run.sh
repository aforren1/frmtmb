#!/usr/bin/env bash
# Release 0.70.0: rerun lane optrap's studies on the merged build
# (rellib-r8), one R process per job, at most P at once, every log in
# dev/rel070-log/studies/optrap/. Environment as dev/optrap-rev3-run.sh,
# with the release tree's stan cache. Each job line:
# <log name> <script> <args...>
#   bash dev/rel070-optrap-run.sh <jobfile> <P>
set -u
cd "$(dirname "$0")/.."
export TMP=/c/Users/adf44/AppData/Local/Temp TEMP=/c/Users/adf44/AppData/Local/Temp
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-release/dev/stan-cache
unset OPTRAP_ABLATE OPTRAP_NOKEEP OPTIMA_NOSEARCH OPTRAP_POLY_ARM OPTRAP_POLY_SD OPTRAP_POLY_YK OPTRAP_REV_RESID OPTRAP_ROUND1
R="${RSCRIPT:-/c/Program Files/R/R-4.6.1/bin/Rscript.exe}"
L=dev/rel070-log/studies/optrap
mkdir -p "$L" "$L/mo"
run1() {
  name="$1"; shift
  echo "START $name $(date +%H:%M:%S) $*" >> "$L/run-ends.txt"
  "$R" "$@" > "$L/$name.txt" 2>&1
  echo "END $name rc=$? $(date +%H:%M:%S)" >> "$L/run-ends.txt"
}
export -f run1
export R L
grep -v '^#' "$1" | grep -v '^$' |
  xargs -P "$2" -L 1 bash -c 'run1 "$@"' _
echo ALLDONE
