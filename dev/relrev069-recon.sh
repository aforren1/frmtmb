#!/usr/bin/env bash
# Reviewer of the 0.69.0 consolidation: rebuild the plain textual merge
# of the four lanes onto 4f5ea39f (order ciharden, surface, setier,
# optima, as dev/rel069-merge.sh does) in a scratch tree, then diff the
# release tree against it. What remains is what the consolidation did
# on top of the lanes: conflict resolutions and its own edits.
# Generated files (Rd, NAMESPACE, the sample copy of helper-brms-suite.R)
# are compared separately by the caller; untracked files no lane has
# are listed as consolidation-new.
set -u
root=/c/Users/adf44/source/r
rel=$root/frmtmb-wt-release
BASE=4f5ea39f
LANES="ciharden surface setier optima"
S=/c/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad
R=$S/relrev069-recon; rm -rf "$R"; mkdir -p "$R"
out=$rel/dev/relrev069-log; mkdir -p "$out"
excl='^(man/.*\.Rd|NAMESPACE|extensions/[^/]*/man/.*\.Rd|extensions/[^/]*/NAMESPACE|extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R|dev/.*-log/.*|dev/stan-cache/.*|dev/relrev069-.*)$'
files=$( (git -C $rel diff --name-only $BASE; git -C $rel ls-files --others --exclude-standard) | grep -vE "$excl" | sort -u)
: > "$out/recon-diff.txt"; : > "$out/recon-new.txt"; : > "$out/recon-conflicts.txt"
for f in $files; do
  mkdir -p "$R/$(dirname "$f")"
  if git -C $rel cat-file -e "$BASE:$f" 2>/dev/null; then
    git -C $rel show "$BASE:$f" | tr -d '\r' > "$R/$f"
    git -C $rel show "$BASE:$f" | tr -d '\r' > "$R/$f.base"
    for l in $LANES; do
      lf=$root/frmtmb-wt-$l/$f
      if [ ! -e "$lf" ]; then
        echo "lane $l deletes $f" >> "$out/recon-conflicts.txt"; continue
      fi
      tr -d '\r' < "$lf" > "$R/$f.theirs"
      cmp -s "$R/$f.theirs" "$R/$f.base" && continue
      git merge-file -L recon -L base -L $l "$R/$f" "$R/$f.base" "$R/$f.theirs" >/dev/null 2>&1 ||
        echo "conflict $l $f" >> "$out/recon-conflicts.txt"
    done
    rm -f "$R/$f.base" "$R/$f.theirs"
  else
    found=""
    for l in $LANES; do
      lf=$root/frmtmb-wt-$l/$f
      if [ -e "$lf" ]; then
        if [ -n "$found" ] && ! cmp -s "$lf" "$R/$f"; then
          echo "added by two lanes, differ: $f ($found $l)" >> "$out/recon-conflicts.txt"
        fi
        cp "$lf" "$R/$f"; found="$found $l"
      fi
    done
    if [ -z "$found" ]; then echo "$f" >> "$out/recon-new.txt"; continue; fi
    # untracked files are copied byte for byte, so compare raw
    if ! cmp -s "$R/$f" "$rel/$f"; then
      echo "=== $f (added by$found)" >> "$out/recon-diff.txt"
      diff -u <(tr -d '\r' < "$R/$f") <(tr -d '\r' < "$rel/$f") | tail -n +3 >> "$out/recon-diff.txt"
    fi
    continue
  fi
  if [ ! -e "$rel/$f" ]; then echo "=== $f DELETED in release" >> "$out/recon-diff.txt"; continue; fi
  tr -d '\r' < "$rel/$f" > "$R/$f.rel"
  if ! cmp -s "$R/$f" "$R/$f.rel"; then
    echo "=== $f" >> "$out/recon-diff.txt"
    diff -u "$R/$f" "$R/$f.rel" | tail -n +3 >> "$out/recon-diff.txt"
  fi
  rm -f "$R/$f.rel"
done
echo "files considered: $(echo "$files" | wc -l)"
echo "differ from recon: $(grep -c '^=== ' "$out/recon-diff.txt")"
grep '^=== ' "$out/recon-diff.txt"
echo "consolidation-new untracked: $(wc -l < "$out/recon-new.txt")"
echo "conflicts in recon:"; cat "$out/recon-conflicts.txt"
