#!/usr/bin/env bash
# Q5 of the codesize lane: size of core R/ at each release commit from
# 0.60.0 to 0.83.0, per file, from read-only git (ls-tree, show).
#   bash dev/codesize-growth.sh
# Writes dev/codesize-out/growth.tsv (version, commit, date, file, bytes,
# lines; bytes are blob sizes, LF line ends) and extracts R/ of 0.60.0
# and 0.83.0 to dev/codesize-out/rel/<version>/R for the line-kind split
# (dev/codesize-growth.R).
set -euo pipefail
export GIT_OPTIONAL_LOCKS=0
wt=/c/Users/adf44/source/r/frmtmb-wt-codesize
out=$wt/dev/codesize-out
cd "$wt"
# version commit pairs: the release commits, found by subject
rels=$(
  git log --format='%h %ad %s' --date=short |
    grep -E '^[0-9a-f]+ [0-9-]+ (Release )?frmtmb 0\.[6-8][0-9]\.[0-9]+[:,]' |
    sed -E 's/^([0-9a-f]+) ([0-9-]+) (Release )?frmtmb (0\.[0-9]+\.[0-9]+).*/\4 \1 \2/' |
    awk '$1 >= "0.60.0"' | sort -t. -k2,2n -k3,3n)
echo "$rels" > "$out/releases.txt"
printf 'version\tcommit\tdate\tfile\tbytes\tlines\n' > "$out/growth.tsv"
while read -r ver com dat; do
  git ls-tree -l "$com" R/ | while read -r mode typ blob size path; do
    n=$(git cat-file blob "$blob" | wc -l)
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$ver" "$com" "$dat" "${path#R/}" \
      "$size" "$n" >> "$out/growth.tsv"
  done
done <<< "$rels"
for v in 0.60.0 0.83.0; do
  com=$(awk -v v="$v" '$1 == v {print $2}' "$out/releases.txt")
  mkdir -p "$out/rel/$v/R"
  git ls-tree --name-only "$com" R/ | while read -r p; do
    git show "$com:$p" > "$out/rel/$v/$p"
  done
  git show "$com:NEWS.md" > "$out/rel/$v/NEWS.md"
done
wc -l < "$out/growth.tsv"
cat "$out/releases.txt"
