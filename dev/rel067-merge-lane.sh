#!/usr/bin/env bash
# Merge one lane's uncommitted changes into the 0.67.0 release worktree
# file by file with git merge-file (base = the lane's HEAD blob, ours =
# the release working file, theirs = the lane's working file). It
# touches no index. Line endings are normalized to LF first, as git
# stores them, so a CRLF checkout against an LF lane file is not a
# whole-file conflict. Generated files (Rd, NAMESPACE, frmtmb.sample's
# helper copy) are left out; untracked files are copied without
# overwriting, except generated Rd and compiled Stan caches.
set -u
lane="$1"
src="/c/Users/adf44/source/r/frmtmb-wt-$lane"
dst="/c/Users/adf44/source/r/frmtmb-wt-release"
S=/c/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad
W="$S/r067-merge-$lane"; rm -rf "$W"; mkdir -p "$W/keep"
base=$(git -C "$src" rev-parse HEAD)
excl='^(man/.*\.Rd|NAMESPACE|extensions/[^/]*/man/.*\.Rd|extensions/[^/]*/NAMESPACE|extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R)$'
git -C "$src" diff --name-status HEAD 2>/dev/null | while read -r st f; do
  echo "$f" | grep -qE "$excl" && continue
  case "$st" in
    M)
      git -C "$src" show "$base:$f" | tr -d '\r' > "$W/base"
      tr -d '\r' < "$src/$f" > "$W/theirs"
      tr -d '\r' < "$dst/$f" > "$W/ours"
      # keep the three inputs of a conflicted file for the hand merge
      k="$W/keep/$(echo "$f" | tr '/' '_')"
      cp "$W/ours" "$k.ours"; cp "$W/base" "$k.base"; cp "$W/theirs" "$k.theirs"
      if git merge-file -L ours -L base -L theirs "$W/ours" "$W/base" "$W/theirs"; then
        r=clean; rm -f "$k.ours" "$k.base" "$k.theirs"
      else
        r="CONFLICT($?)"
      fi
      cp "$W/ours" "$dst/$f"
      echo "   $r $f";;
    A)
      if [ -e "$dst/$f" ]; then echo "   CLASH added: $f"; else
        mkdir -p "$dst/$(dirname "$f")"; cp -p "$src/$f" "$dst/$f"; echo "   added $f"; fi;;
    *) echo "   UNHANDLED $st $f";;
  esac
done
lst="$S/r067-lane-$lane.untracked"; : > "$lst"; clash=0; skip=0
while IFS= read -r f; do
  case "$f" in
    man/*.Rd|extensions/*/man/*.Rd|*stan-cache*/*) skip=$((skip+1)); continue;;
  esac
  if [ -e "$dst/$f" ]; then echo "   CLASH untracked: $f"; clash=$((clash+1)); continue; fi
  printf '%s\n' "$f" >> "$lst"
done < <(git -C "$src" ls-files --others --exclude-standard)
( cd "$src" && tar -cf - -T "$lst" ) | ( cd "$dst" && tar -xpf - )
echo "   untracked copied: $(wc -l < "$lst"), skipped: $skip, clashes: $clash"
