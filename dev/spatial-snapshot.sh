#!/usr/bin/env bash
# Write the lane's whole change (tracked diff plus every untracked file
# outside dev/) as one patch, dev/spatial-out/<tag>.diff, read-only on
# git: nothing is staged.
set -e
cd "$(dirname "$0")/.."
out="dev/spatial-out/$1.diff"
mkdir -p dev/spatial-out
git -c core.safecrlf=false diff --binary > "$out" 2>/dev/null
git ls-files --others --exclude-standard | grep -v '^dev/' | while read -r f; do
  git -c core.safecrlf=false diff --no-index --binary /dev/null "$f" >> "$out" 2>/dev/null || true
done
echo "$out: $(grep -c '^diff --git' "$out") files"
# and a copy of every changed or new file outside dev/, for the
# stage-to-stage patches (dev/spatial-stage-diff.sh)
tree="dev/spatial-out/tree-$1"
rm -rf "$tree"; mkdir -p "$tree"
{ git diff --name-only; git ls-files --others --exclude-standard; } |
  grep -v '^dev/' | sort -u | while read -r f; do
  [ -f "$f" ] || continue
  mkdir -p "$tree/$(dirname "$f")"; cp "$f" "$tree/$f"
done
echo "$tree: $(find "$tree" -type f | wc -l) files"
