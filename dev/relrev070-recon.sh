#!/usr/bin/env bash
# Reviewer of the 0.70.0 consolidation: rebuild the plain textual merge of
# the three lanes onto 2884ccbd, independent of dev/rel070-merge.sh, so the
# consolidation's own edits can be told apart from what the lanes wrote.
set -u
R=C:/Users/adf44/source/r
REL=$R/frmtmb-wt-release
OUT=$REL/dev/relrev070-out/recon
rm -rf "$OUT"; mkdir -p "$OUT"
LOG=$REL/dev/relrev070-log/recon.txt
mkdir -p "$(dirname "$LOG")"; : > "$LOG"
files=$( (for l in moparity optrap seflags; do
  git -C $R/frmtmb-wt-$l diff --name-only 48c1e6c9 2>/dev/null; done) |
  sort -u | grep -v '^man/' | grep -v '/man/' | grep -v '^NAMESPACE$')
for f in $files; do
  mkdir -p "$OUT/$(dirname "$f")"
  tmp=$OUT/.tmp; mkdir -p "$tmp"
  git -C "$REL" show "2884ccbd:$f" | tr -d '\r' > "$OUT/$f"
  git -C "$REL" show "48c1e6c9:$f" | tr -d '\r' > "$tmp/base"
  for l in moparity optrap seflags; do
    if git -C $R/frmtmb-wt-$l diff --quiet 48c1e6c9 -- "$f" 2>/dev/null; then
      continue
    fi
    tr -d '\r' < "$R/frmtmb-wt-$l/$f" > "$tmp/theirs"
    git merge-file -L release -L base -L "$l" "$OUT/$f" "$tmp/base" \
      "$tmp/theirs"
    echo "merged $l $f conflicts=$?" >> "$LOG"
  done
  if [ -f "$REL/$f" ]; then
    if diff -q <(tr -d '\r' < "$REL/$f") "$OUT/$f" > /dev/null; then
      echo "SAME $f" >> "$LOG"
    else
      n=$(diff <(tr -d '\r' < "$OUT/$f") <(tr -d '\r' < "$REL/$f") | grep -c '^[<>]')
      echo "DIFF $f lines=$n" >> "$LOG"
    fi
  else
    echo "MISSING-IN-RELEASE $f" >> "$LOG"
  fi
done
rm -rf "$OUT/.tmp"
# Files the release changed against 2884ccbd that no lane changed.
git -C "$REL" diff --name-only 2884ccbd | while read -r f; do
  if ! echo "$files" | grep -qx "$f"; then echo "RELONLY $f" >> "$LOG"; fi
done
cat "$LOG"
