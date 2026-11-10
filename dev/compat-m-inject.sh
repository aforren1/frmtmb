#!/usr/bin/env bash
# Run the fixed dev/equiv-measure.R in a copy of the instrument whose
# dev/equiv-live-row.R is dev/compat-m-inject-row.R, once per failure
# mode, one pair, one job at a time. The worktree's files are not
# touched.
#   bash dev/compat-m-inject.sh <id> <mode> ...
set -u
id=$1; shift
W="$(cd "$(dirname "$0")/.." && pwd)"
T="$W/dev/compat-m-inject"
rm -rf "$T"; mkdir -p "$T/dev/equiv-data"
for f in equiv-measure.R equiv-lib.R equiv-specs.R equiv-flow.R \
         equiv-dump-compat.R; do cp "$W/dev/$f" "$T/dev/$f"; done
cp "$W/dev/equiv-data/sample.tsv" "$T/dev/equiv-data/"
mkdir -p "$T/tests/testthat"
for f in helper-brms.R helper-stan.R helper-brms-methods.R; do
  cp "$W/tests/testthat/$f" "$T/tests/testthat/$f"; done
cp "$W/dev/equiv-live-row.R" "$T/dev/equiv-live-row-real.R"
cp "$W/dev/compat-m-inject-row.R" "$T/dev/equiv-live-row.R"
export FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-compat-lib;C:/Users/adf44/source/r/rellib-r19"
export FRMTMB_STAN_CACHE="$W/dev/stan-cache"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd "$T"
for m in "$@"; do
  COMPAT_INJECT=$m "$RS" dev/equiv-measure.R "$id" "out-$m" --force \
    > "log-$m.txt" 2>&1
  echo "== mode $m: exit status $?"
  cat "log-$m.txt"
done
