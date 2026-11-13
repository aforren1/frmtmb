#!/usr/bin/env bash
# Re-merge lane boundfam's consolidation fixes (its findings, section
# "Consolidation fixes", 2026-11-12) into the release tree. Every path
# boundfam takes now is compared with the release tree, line endings
# aside. A path only boundfam takes is copied when it differs; a path
# other lanes take too is merged three-way against e8abed0a (the release
# tree already holds boundfam's earlier hunks, which merge as identical
# changes). NEWS.md and dev/test-backlog.md are rebuilt by
# dev/rel082-news.R and dev/rel082-remerge-one.sh + dev/rel082-backlog.R.
# Written back with the release checkout's line endings.
#
#   bash dev/rel082-remerge-boundfam.sh
set -euo pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
REL=$R/frmtmb-wt-rel082
B=$R/frmtmb-wt-boundfam
BASE=e8abed0a
LANES="quadseam ppdraw wcheck priortab nlstop gpstop ordstart silent arspat names sampfront"
OUT=$REL/dev/rel082-log/remerge-boundfam
rm -rf "$OUT"; mkdir -p "$OUT/tmp"
lst() { (cd "$1" && git status --porcelain=v1 -z --untracked-files=all --no-renames) |
          tr '\0' '\n' | cut -c4- | sort -u; }
lst "$B" > "$OUT/boundfam.txt"
: > "$OUT/others.txt"
for l in $LANES; do lst "$R/frmtmb-wt-$l" >> "$OUT/others.txt"; done
sort -u -o "$OUT/others.txt" "$OUT/others.txt"
SKIP='^(NEWS\.md|dev/test-backlog\.md)$'
crlf() { sed 's/$/\r/' | perl -0pe 's/\r\z//'; }
n=0
while IFS= read -r p; do
  [ -n "$p" ] || continue
  [[ "$p" =~ $SKIP ]] && continue
  if [ -e "$REL/$p" ] && cmp -s <(tr -d '\r' < "$B/$p") <(tr -d '\r' < "$REL/$p"); then
    continue
  fi
  if [ -e "$REL/$p" ] && [ "$(tr -dc '\r' < "$REL/$p" | head -c 1 | wc -c)" -gt 0 ]; then
    eol=crlf
  else
    eol=$( [ -e "$REL/$p" ] && echo lf || echo new )
  fi
  if grep -qxF -- "$p" "$OUT/others.txt" && git -C "$REL" cat-file -e "$BASE:$p" 2>/dev/null; then
    # The release copy holds boundfam's earlier hunks, so a three-way
    # merge against e8abed0a conflicts where boundfam rewrote its own
    # lines. The file is rebuilt from the base and every lane instead, in
    # the merge's order, which holds only while the consolidation has not
    # edited it after the merge.
    case " R/fit.R R/predict-brms.R R/autoscale.R R/gp-start.R " in
      *" $p "*) echo "STOP: $p was edited by the consolidation"; exit 1 ;;
    esac
    bash "$REL/dev/rel082-remerge-one.sh" "$p" > "$OUT/tmp/remerge.txt" 2>&1 ||
      { echo "CONFLICT $p"; cat "$OUT/tmp/remerge.txt"; exit 1; }
    if grep -qE '^(<<<<<<<|>>>>>>>) ' "$REL/$p"; then echo "MARKERS in $p"; exit 1; fi
    echo "REBUILT $p ($(grep -c '^merged' "$OUT/tmp/remerge.txt") lanes)" | tee -a "$OUT/done.txt"
  else
    mkdir -p "$(dirname "$REL/$p")"
    if [ "$eol" = crlf ]; then tr -d '\r' < "$B/$p" | crlf > "$REL/$p"; else cp -p "$B/$p" "$REL/$p"; fi
    echo "COPIED $p ($eol)" | tee -a "$OUT/done.txt"
  fi
  n=$((n + 1))
done < "$OUT/boundfam.txt"
echo "REMERGE DONE: $n paths; boundfam takes $(grep -c . "$OUT/boundfam.txt")"
