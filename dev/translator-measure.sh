#!/usr/bin/env bash
# Measure a list of survey pairs into one outdir, at most `jobs` at once
# (the round allows three R processes per lane):
#
#   bash dev/translator-measure.sh <outdir> <jobs> <id> [<id> ...]
#
# from the worktree root. The library list comes from FRMTMB_EQUIV_LIB,
# as dev/equiv-measure.R reads it. Every pair is forced, so a rerun into
# the same outdir remeasures. Never writes the record.
set -u
out=$1
jobs=$2
shift 2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="$(cd dev/stan-cache && pwd -W)"
export NOT_CRAN=true
mkdir -p "$out/logs"
start=$(date +%s)
printf '%s\n' "$@" | xargs -P "$jobs" -I{} bash -c \
  "timeout 2400 \"$R\" dev/equiv-measure.R {} \"$out\" --force \
     > \"$out/logs/{}.log\" 2>&1 ||
   echo \"FAILED {} exit \$?\" >> \"$out/logs/{}.log\""
end=$(date +%s)
echo "measured $# pairs in $((end - start)) s"
