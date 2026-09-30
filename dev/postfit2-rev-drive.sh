#!/usr/bin/env bash
# Reviewer driver: run the jobs listed in a file (one per line:
# "<tag> <pkg> <file> <arm> [flags]"), P at a time, one R process per
# job, log to dev/postfit2-rev-log/<tag>.txt.
#   /usr/bin/bash dev/postfit2-rev-drive.sh <jobfile> [P]
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-postfit2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$WT/dev/postfit2-rev-log
mkdir -p "$OUT"
P=${2:-16}
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export TMP="C:/Users/adf44/AppData/Local/Temp" TEMP="C:/Users/adf44/AppData/Local/Temp"
grep -v '^#' "$1" | grep -v '^$' | xargs -P "$P" -L 1 /usr/bin/bash -c \
  'tag="$0"; "'"$R"'" '"$WT"'/dev/postfit2-rev-runtest.R "$@" > '"$OUT"'/"$tag".txt 2>&1'
echo "jobs: $(grep -v '^#' "$1" | grep -vc '^$')"
