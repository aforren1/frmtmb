#!/usr/bin/env bash
# Lane optrap, punch round 2: run a job list, one R process each, at
# most P at once, with every log in dev/optrap-log/p2/. Each line:
# <log name> <script> <args...>; a reviewer's script (dev/optrap-rev*)
# runs through dev/optrap-p2-redirect.R so its outputs land in
# dev/optrap-log/rev/ and dev/optrap-log/rev2/.
#   bash dev/optrap-p2-run.sh <jobfile> <P>
set -u
cd "$(dirname "$0")/.."
export TMP=/c/Users/adf44/AppData/Local/Temp TEMP=/c/Users/adf44/AppData/Local/Temp
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-optrap/dev/stan-cache
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/optrap-log/p2 dev/optrap-log/p2/mo
run1() {
  name="$1"; script="$2"; shift 2
  case "$script" in
    dev/optrap-rev*) "$R" dev/optrap-p2-redirect.R "$script" "$@" ;;
    *) "$R" "$script" "$@" ;;
  esac > "dev/optrap-log/p2/$name.txt" 2>&1
  echo "END $name rc=$?" >> dev/optrap-log/p2/run-ends.txt
}
export -f run1
export R
grep -v '^#' "$1" | grep -v '^$' |
  xargs -P "$2" -L 1 bash -c 'run1 "$@"' _
echo ALLDONE
