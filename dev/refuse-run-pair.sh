#!/usr/bin/env bash
# Measure sampled pairs on this lane's build, into dev/refuse-pair-out,
# with --force, so a pair can be measured again while a fix moves. The
# rerun against the baseline is dev/equiv-run.sh with dev/refuse-equiv-out.
#   bash dev/refuse-run-pair.sh p043 p070 ...
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="$(cd dev/stan-cache && pwd -W)"
export NOT_CRAN=true
export FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-refuse-lib;C:/Users/adf44/source/r/rellib-r16"
out=dev/refuse-pair-out
mkdir -p "$out/logs"
printf '%s\n' "$@" | xargs -P 8 -I{} bash -c \
  "timeout 2400 \"$R\" dev/equiv-measure.R {} \"$out\" --force \
     > \"$out/logs/{}.log\" 2>&1 || echo \"FAILED {} exit \$?\" >> \"$out/logs/{}.log\""
for id in "$@"; do
  echo "== $id"
  "$R" -e "r <- read.delim('$out/$id.tsv'); keep <- intersect(c('verdict','level','l1_const','l1_grad','l1_extra_const','l1_error','frm_error','brms_error','l2','reason'), names(r)); for (k in keep) cat(k, ':', substr(as.character(r[[k]]), 1, 300), '\n')" 2>&1
done
