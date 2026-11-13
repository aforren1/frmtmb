#!/usr/bin/env bash
# Run one dev/quadseam-*.R script against a library tag, logging to
# dev/quadseam-log/<name>-<tag>.txt. Usage:
#   bash dev/quadseam-run.sh <script-name> <base|lib> [extra args]
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam
LIB=C:/Users/adf44/source/r/wt-quadseam-lib
name=$1; tag=$2; shift 2
suf=${QS_SUF:-}
arg=$tag
[ "$tag" = lib ] && arg=$LIB
mkdir -p "$WT/dev/quadseam-log"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
cd "$WT"
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "dev/quadseam-$name.R" "$arg" "$@" \
  > "dev/quadseam-log/$name-$tag$suf.txt" 2>&1
echo "exit $? $name $tag"
