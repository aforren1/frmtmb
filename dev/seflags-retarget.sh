#!/usr/bin/env bash
# Lane seflags: copy the SE attack-set scripts of lanes nanse and setier
# (and the 0.69.0 review's) into dev/seflags-out/scripts/, pointed at
# this lane's libraries, so that the originals stay as their lanes left
# them. "lane" arms read wt-seflags-lib over rellib-r7; "base" arms and
# every rellib-r6 read rellib-r7 (frmtmb 0.69.0); helper paths into the
# release worktree read this worktree's tests; a table a script wrote
# into dev/setier-rev2-log/ goes to $SEFLAGS_OUT (default
# dev/seflags-log).
#
#   bash dev/seflags-retarget.sh <script.R> ...
set -eu
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
OUT="$W/dev/seflags-out/scripts"
mkdir -p "$OUT"
for s in "$@"; do
  b=$(basename "$s")
  sed -e 's#C:/Users/adf44/source/r/wt-setier-lib#C:/Users/adf44/source/r/wt-seflags-lib#g' \
      -e 's#C:/Users/adf44/source/r/wt-optima-lib#C:/Users/adf44/source/r/wt-seflags-lib#g' \
      -e 's#C:/Users/adf44/source/r/rellib-r6#C:/Users/adf44/source/r/rellib-r7#g' \
      -e 's#C:/Users/adf44/source/r/frmtmb-wt-release#C:/Users/adf44/source/r/frmtmb-wt-seflags#g' \
      -e 's#C:/Users/adf44/source/r/frmtmb-wt-setier#C:/Users/adf44/source/r/frmtmb-wt-seflags#g' \
      -e 's#"dev/setier-rev2-log"#Sys.getenv("SEFLAGS_OUT", "dev/seflags-log")#g' \
      -e 's#"dev/setier-rev2-log/#paste0(Sys.getenv("SEFLAGS_OUT", "dev/seflags-log"), "/"), "#g' \
      "$W/$s" > "$OUT/$b"
done
