#!/usr/bin/env bash
# Reviewer: rebuild the textual merge of the four lanes onto 9e902909
# (order vigport, fixes, gpby, ordmix, as the integration did) in a
# scratch tree, then diff the release tree against it. What remains is
# what the consolidation did by hand: conflict resolutions and edits.
set -u
root=/c/Users/adf44/source/r
rel=$root/frmtmb-wt-release
S=/c/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad
R=$S/relrev-recon; rm -rf "$R"; mkdir -p "$R"
out=$rel/dev/relrev-log; mkdir -p "$out"
excl='^(man/.*\.Rd|NAMESPACE|extensions/[^/]*/man/.*\.Rd|extensions/[^/]*/NAMESPACE|extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R|dev/.*)$'
files=$( (git -C $rel diff --name-only 9e902909; git -C $rel ls-files --others --exclude-standard) | grep -vE "$excl" | sort -u)
: > "$out/recon-diff.txt"
for f in $files; do
  mkdir -p "$R/$(dirname "$f")"
  if git -C $rel cat-file -e "9e902909:$f" 2>/dev/null; then
    git -C $rel show "9e902909:$f" | tr -d '\r' > "$R/$f"
    git -C $rel show "9e902909:$f" | tr -d '\r' > "$R/$f.base"
    for l in vigport fixes gpby ordmix; do
      lf=$root/frmtmb-wt-$l/$f
      [ -e "$lf" ] || continue
      tr -d '\r' < "$lf" > "$R/$f.theirs"
      cmp -s "$R/$f.theirs" "$R/$f.base" && continue
      git merge-file -L recon -L base -L $l "$R/$f" "$R/$f.base" "$R/$f.theirs" >/dev/null 2>&1
    done
    rm -f "$R/$f.base" "$R/$f.theirs"
  else
    # added file: take the one lane that has it (or none)
    found=""
    for l in vigport fixes gpby ordmix; do
      lf=$root/frmtmb-wt-$l/$f
      if [ -e "$lf" ]; then tr -d '\r' < "$lf" > "$R/$f"; found="$found $l"; fi
    done
    [ -z "$found" ] && : > "$R/$f"
  fi
  tr -d '\r' < "$rel/$f" > "$R/$f.rel"
  if ! cmp -s "$R/$f" "$R/$f.rel"; then
    echo "=== $f" >> "$out/recon-diff.txt"
    diff -u "$R/$f" "$R/$f.rel" | tail -n +3 >> "$out/recon-diff.txt"
  fi
  rm -f "$R/$f.rel"
done
grep -c '^=== ' "$out/recon-diff.txt"
grep '^=== ' "$out/recon-diff.txt"
