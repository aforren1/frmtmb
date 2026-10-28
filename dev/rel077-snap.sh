#!/usr/bin/env bash
# md5 of every tracked and untracked (not ignored) file of the release
# tree, outside dev/ into <out> and under dev/ into <out>.dev, to show
# that a rerun of dev/rel077-merge.sh changed nothing it should not.
#   bash dev/rel077-snap.sh <out>
cd /c/Users/adf44/source/r/frmtmb-wt-rel077 || exit 1
{ git ls-files; git ls-files --others --exclude-standard; } | sort -u |
  grep -v '^dev/' | while read -r f; do [ -f "$f" ] && md5sum "$f"; done > "$1"
{ git ls-files dev; git ls-files --others --exclude-standard dev; } |
  sort -u | while read -r f; do [ -f "$f" ] && md5sum "$f"; done > "$1.dev"
