#!/usr/bin/env bash
# Reviewer of lane optrap, final check: run a list of jobs, one R process each, at
# most P at once. Each line of the job file: <log name> <script> <args...>
#   bash dev/optrap-rev-run.sh <jobfile> <P>
set -u
cd "$(dirname "$0")/.."
export TMP=/c/Users/adf44/AppData/Local/Temp TEMP=/c/Users/adf44/AppData/Local/Temp
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-optrap/dev/stan-cache
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/optrap-rev3-log
run1() {
  name="$1"; shift
  "$R" "$@" > "dev/optrap-rev3-log/$name.txt" 2>&1
  echo "END $name rc=$?" >> dev/optrap-rev3-log/run-ends.txt
}
export -f run1
export R
grep -v '^#' "$1" | grep -v '^$' |
  xargs -P "$2" -L 1 bash -c 'run1 "$@"' _
echo ALLDONE
