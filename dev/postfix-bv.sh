#!/usr/bin/env bash
# Lane postfix: the hand-translated brms vignettes this lane's items
# reach, against one build, one R process per script.
#
#   bash dev/postfix-bv.sh <arm> [vignettes...]
#
# arm "before" reads rellib-r8 alone, "after" the lane library in front
# of it. One CSV per script lands in dev/postfix-bv-out/<arm>/ (BV_OUT),
# with one log per script.
set -u
arm="$1"; shift
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
here="$ROOT/dev/brms-vignettes"
out="$ROOT/dev/postfix-bv-out/$arm"
mkdir -p "$out"
rscript="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
if [ "$arm" = "before" ]; then
  lib="C:/Users/adf44/source/r/rellib-r8"
else
  lib="C:/Users/adf44/source/r/wt-postfix-lib;C:/Users/adf44/source/r/rellib-r8"
fi
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export BV_LIB="$lib" BV_OUT="$out" BV_DIR="$here"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
for nm in "$@"; do
  "$rscript" "$here/$nm.R" > "$out/bv-$nm.log" 2>&1 &
done
wait
for nm in "$@"; do
  f="$out/bv-$nm.log"
  printf '%-30s %s ok, %s err, %s\n' "$nm" "$(grep -c '^\[ok' "$f")" \
    "$(grep -c '^\[err' "$f")" "$(grep -m1 '^### done' "$f" || echo NO DONE)"
done
