#!/usr/bin/env bash
# Reruns the review's scripts on the current build, to the lane's own logs.
cd /c/Users/adf44/source/r/frmtmb-wt-nlfix || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
: > dev/nlfix-rc.progress
run() {
  arm=$1; s=$2; shift 2
  FRM_ARM=$arm "$R" "dev/review-nlfix-$s.R" "$@" > "dev/nlfix-rc-$s$(printf '%s' "$*")-$arm.log" 2>&1
  echo "$(date +%H:%M:%S) $s $* $arm $?" >> dev/nlfix-rc.progress
}
run after c1b
run after c1
run after c3b
run after c4
run before c4
run after c5 bbgrid
run after c3
run before c3b
