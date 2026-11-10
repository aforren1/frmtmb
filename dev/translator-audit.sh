#!/usr/bin/env bash
# Run dev/translator-audit.R over a list of pairs, at most `jobs` at once:
#
#   bash dev/translator-audit.sh <outdir> <jobs> <id> [<id> ...]
#
# from the worktree root, with FRMTMB_EQUIV_LIB set as for
# dev/equiv-measure.R.
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
  "timeout 1200 \"$R\" dev/translator-audit.R {} \"$out\" \
     > \"$out/logs/{}.log\" 2>&1 ||
   echo \"FAILED {} exit \$?\" >> \"$out/logs/{}.log\""
end=$(date +%s)
echo "audited $# pairs in $((end - start)) s"
