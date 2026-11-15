#!/usr/bin/env bash
# Review of lane quadseam2: run dev/reviews/quadseam2-review-<name>.R on
# one library, pinned to one P-core mask (rtmb-pitfalls.md item 21).
#   bash quadseam2-review-run.sh <name> <base|lib> <mask> [suffix] [args]
# Log: dev/reviews/quadseam2-review-log/<name>-<tag><suffix>.txt
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam2
LIB=C:/Users/adf44/source/r/wt-quadseam2-lib
name=$1; tag=$2; mask=$3; suf=${4:-}; shift 3; [ $# -gt 0 ] && shift
arg=$tag
# the R script maps "lib" to the lane library (Git Bash mangles paths)
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_PARALLEL_PORT=$((11000 + RANDOM % 5000))
cd "$WT"
cmd //c "dev\\reviews\\quadseam2-review-aff.cmd $mask dev/reviews/quadseam2-review-$name.R $arg $*" \
  > "dev/reviews/quadseam2-review-log/$name-$tag$suf.txt" 2>&1
echo "exit $? $name $tag$suf"
