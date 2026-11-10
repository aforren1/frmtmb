#!/usr/bin/env bash
# Measure a list of pairs with the instrument as it stood at e0cce906
# (dev/compat-equiv-measure-orig.R, a copy) and with the edited
# dev/equiv-measure.R, on one library, COMPAT_JOBS R processes at once
# (default 3; the edited arm starts a child R for its live row, so use 1
# to stay within three).
#   bash dev/compat-noreg-run.sh <lib;lib> <outprefix> <id> ...
set -u
export FRMTMB_EQUIV_LIB=$1
pre=$2
shift 2
cd "$(dirname "$0")/.."
export FRMTMB_STAN_CACHE="$PWD/dev/stan-cache"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p "$pre-orig-out" "$pre-edit-out"
n=0
for id in "$@"; do
  for arm in orig edit; do
    if [ "$arm" = orig ]; then s=dev/compat-equiv-measure-orig.R; else
      s=dev/equiv-measure.R; fi
    ( "$RS" "$s" "$id" "$pre-$arm-out" --force > "$pre-$arm-out/$id.log" 2>&1 ) &
    n=$((n + 1))
    if [ $((n % ${COMPAT_JOBS:-3})) -eq 0 ]; then wait; fi
  done
done
wait
echo "done $n runs"
