#!/usr/bin/env bash
# Lane sparsese: run each line of a jobs file (arguments to Rscript, then
# `>` and a log path relative to the worktree) at most two at a time.
#   bash dev/sparsese-jobs.sh <jobs-file>
set -u
ROOT="C:/Users/adf44/source/r/frmtmb-wt-sparsese"
cd "$ROOT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/sparsese-out
while IFS= read -r line; do
  [ -z "$line" ] && continue
  args="${line%% > *}"
  log="${line##* > }"
  while [ "$(jobs -rp | wc -l)" -ge 2 ]; do wait -n; done
  echo "start $(date +%T) $args"
  # shellcheck disable=SC2086
  ( "$R" $args > "$log" 2>&1; echo "exit $? $(date +%T) $log" ) &
done < "$1"
wait
echo "all done $(date +%T)"
