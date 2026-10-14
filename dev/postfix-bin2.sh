#!/usr/bin/env bash
# Lane postfix: lane remeasure's bin-2 port of brms's suite, run from
# this worktree (its helpers, staged into dev/brmsport-bin2-run/) on
# rellib-r8 and on the lane library, one file at a time because the
# runner stages into one directory per package. Records land in
# dev/postfix-bin2-log/<pkg>-<file>-<arm>.tsv. Read-only on remeasure.
#   bash dev/postfix-bin2.sh <arm> <pkg> <file stem>...
set -u
arm="$1"; pkg="$2"; shift 2
RM=/c/Users/adf44/source/r/frmtmb-wt-remeasure/dev
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export FRMTMB_STAN_CACHE="$PWD/dev/stan-cache"
if [ "$arm" = before ]; then
  export FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r8
else
  export FRMTMB_PORT_LIB=C:/Users/adf44/source/r/wt-postfix-lib
fi
mkdir -p dev/postfix-bin2-log
for s in "$@"; do
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$RM/brmsport-bin2-run.R" \
    "$pkg" "$RM/brmsport-bin2/$pkg/test-brms-suite-bin2-$s.R" \
    "dev/postfix-bin2-log/$pkg-$s-$arm.tsv" \
    > "dev/postfix-bin2-log/$pkg-$s-$arm.txt" 2>&1
  tail -1 "dev/postfix-bin2-log/$pkg-$s-$arm.txt"
done
