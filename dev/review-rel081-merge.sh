#!/usr/bin/env bash
# Release review of 0.81.0: is every lane path in the release tree, and
# where the release file differs from the lane's, which lanes touched it.
# Read-only on every tree. Output: dev/review-rel081-merge.txt
set -u
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
REL=$R/frmtmb-wt-rel081
OUT=$REL/dev/review-rel081-merge.txt
EXCL='^dev/(review-[^/]*-out[^/]*|[^/]*-equiv-out[^/]*|review-[^/]*-merge|compat-m-race[^/]*)/'
: > "$OUT"
declare -A owners
for l in translator compat code8 hollow support; do
  (cd "$R/frmtmb-wt-$l" &&
     git status --porcelain=v1 -z --untracked-files=all --no-renames |
     tr '\0' '\n' | cut -c4- | sort -u) > "$REL/dev/review-rel081-paths-$l.lst"
  while read -r p; do
    owners[$p]="${owners[$p]:-} $l"
  done < "$REL/dev/review-rel081-paths-$l.lst"
done
for p in "${!owners[@]}"; do
  if printf '%s\n' "$p" | grep -Eq "$EXCL"; then continue; fi
  set -- ${owners[$p]}
  n=$#
  for l in "$@"; do
    lf=$R/frmtmb-wt-$l/$p
    rf=$REL/$p
    if [ ! -e "$lf" ]; then
      if [ -e "$rf" ]; then echo "DELETED-IN-LANE-PRESENT-IN-REL $l $p" >> "$OUT"; fi
      continue
    fi
    if [ ! -e "$rf" ]; then echo "MISSING $l $p" >> "$OUT"; continue; fi
    if ! cmp -s <(tr -d '\r' < "$lf") <(tr -d '\r' < "$rf"); then
      echo "DIFF n=$n [${owners[$p]# }] $l $p" >> "$OUT"
    fi
  done
done
sort -o "$OUT" "$OUT"
echo "paths: ${#owners[@]}" >> "$OUT"
