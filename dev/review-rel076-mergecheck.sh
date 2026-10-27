#!/usr/bin/env bash
# Review: an independent replay of the merge, read-only on every git
# tree (status with --no-optional-locks, show, merge-file -p to stdout).
# For each path any lane changed (tracked or untracked, generated man/
# and NAMESPACE included), build the expected release file: the lane's
# copy if one lane changed it, the three-way merge of the lane copies
# over 70476dc3 if two or three did. Compare with the release tree,
# LF-normalized. Every difference is printed with its diff, for the
# reviewer to match to a recorded resolution or release fix.
#   bash dev/review-rel076-mergecheck.sh > dev/review-rel076-out/mergecheck.txt
R=/c/Users/adf44/source/r
REL=$R/frmtmb-wt-rel076
W=$REL/dev/review-rel076-out/mc
BASE=70476dc3
rm -rf "${W:?}"; mkdir -p "$W"
lanes="mocode acprior loglik"
for l in $lanes; do
  git -C "$R/frmtmb-wt-$l" --no-optional-locks status --porcelain -uall |
    sed 's/^...//' | grep -v '^dev/review-' | sort > "$W/$l.paths"
  echo "lane $l: $(wc -l < "$W/$l.paths") paths (dev/review-* excluded)"
done
sort -u "$W"/*.paths > "$W/all.paths"
echo "union: $(wc -l < "$W/all.paths") paths"
lf() { tr -d '\r' < "$1"; }
same=0; differ=0; missing=0
while read -r p; do
  ls_=""; for l in $lanes; do grep -qxF "$p" "$W/$l.paths" && ls_="$ls_ $l"; done
  set -- $ls_
  k=$(echo "$p" | tr '/' '_')
  if [ ! -f "$REL/$p" ]; then
    # a lane deleted it, or the release tree did not take it
    if git -C "$REL" cat-file -e "$BASE:$p" 2>/dev/null; then
      echo "ABSENT in release (tracked at base): $p [$ls_ ]"
    else
      echo "ABSENT in release (lane-new): $p [$ls_ ]"
    fi
    missing=$((missing + 1)); continue
  fi
  if [ $# -eq 1 ]; then
    if [ -f "$R/frmtmb-wt-$1/$p" ]; then lf "$R/frmtmb-wt-$1/$p" > "$W/$k.exp"
    else echo "lane $1 deleted $p"; continue; fi
  else
    if git -C "$REL" cat-file -e "$BASE:$p" 2>/dev/null; then
      git -C "$REL" show "$BASE:$p" | tr -d '\r' > "$W/$k.base"
    else : > "$W/$k.base"; fi
    lf "$R/frmtmb-wt-$1/$p" > "$W/$k.exp"
    shift
    for l in "$@"; do
      lf "$R/frmtmb-wt-$l/$p" > "$W/$k.theirs"
      git merge-file -p -L ours -L base -L "$l" "$W/$k.exp" "$W/$k.base" \
        "$W/$k.theirs" > "$W/$k.tmp"
      mv "$W/$k.tmp" "$W/$k.exp"
    done
  fi
  lf "$REL/$p" > "$W/$k.rel"
  if cmp -s "$W/$k.exp" "$W/$k.rel"; then same=$((same + 1))
  else
    differ=$((differ + 1))
    echo "=== DIFFERS: $p [$ls_ ] (< expected from lanes, > release)"
    diff "$W/$k.exp" "$W/$k.rel" | head -60
  fi
done < "$W/all.paths"
echo "SUMMARY same=$same differ=$differ absent=$missing"
# release-tree changes that no lane made (the consolidation's own)
git -C "$REL" --no-optional-locks status --porcelain -uall | sed 's/^...//' |
  grep -v '^dev/review-rel076' | sort > "$W/rel.paths"
echo "release-only paths: $(comm -23 "$W/rel.paths" "$W/all.paths" | wc -l)"
comm -23 "$W/rel.paths" "$W/all.paths" | sed 's/^/  /'
