#!/usr/bin/env bash
# Run one dev/quadseam2-*.R script against a library tag, logging to
# dev/quadseam2-log/<name>-<tag><QS_SUF>.txt. Usage:
#   bash dev/quadseam2-run.sh <script-name> <base|lib> [extra args]
# "base" is rellib-r21 alone; "lib" puts wt-quadseam2-lib in front of it.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam2
LIB=C:/Users/adf44/source/r/wt-quadseam2-lib
name=$1; tag=$2; shift 2
suf=${QS_SUF:-}
arg=$tag
[ "$tag" = lib ] && arg=$LIB
mkdir -p "$WT/dev/quadseam2-log"
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_PARALLEL_PORT=${R_PARALLEL_PORT:-$((11000 + RANDOM % 5000))}
cd "$WT"
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "dev/quadseam2-$name.R" "$arg" "$@" \
  > "dev/quadseam2-log/$name-$tag$suf.txt" 2>&1
echo "exit $? $name $tag"
