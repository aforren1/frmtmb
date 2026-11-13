#!/usr/bin/env bash
# Redo the merge of one file of the release tree from the base and the
# twelve lanes, as dev/rel082-merge.sh does it (union for NEWS files and
# the backlog, three-way otherwise), and write it back CRLF. Used once,
# after the first run of dev/rel082-backlog.R dropped the backlog's blank
# lines (strsplit() of ""), to rebuild the merged backlog before its rerun.
#
#   bash dev/rel082-remerge-one.sh dev/test-backlog.md
set -euo pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
REL=$R/frmtmb-wt-rel082
BASE=e8abed0a
ORDER="boundfam quadseam ppdraw wcheck priortab nlstop gpstop ordstart silent arspat names sampfront"
p=${1:?path}
T=$(mktemp -d)
git -C "$REL" show "$BASE:$p" > "$T/base"
cp "$T/base" "$T/m"
for l in $ORDER; do
  L=$R/frmtmb-wt-$l
  if [ -n "$(git -C "$L" status --porcelain -- "$p")" ]; then
    tr -d '\r' < "$L/$p" > "$T/theirs"
    case "$p" in
      NEWS.md|extensions/*/NEWS.md|dev/test-backlog.md|dev/upstream-bugs.md)
        git merge-file --union -L rel -L base -L "$l" "$T/m" "$T/base" "$T/theirs" ;;
      *) git merge-file -L rel -L base -L "$l" "$T/m" "$T/base" "$T/theirs" ;;
    esac
    echo "merged $l"
  fi
done
sed 's/$/\r/' "$T/m" | perl -0pe 's/\r\z//' > "$REL/$p"
rm -rf "$T"
echo "REMERGED $p"
