#!/usr/bin/env bash
# Release review of 0.81.0: for each file two or more lanes changed, is
# every line a lane added present in the release tree, and every line a
# lane removed absent (unless the base holds it more than once)? Then a
# fresh three-way merge of the lanes' files, LF-normalized, against the
# release file. Read-only. Output dev/review-rel081-hunks.txt
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
REL=$R/frmtmb-wt-rel081
T=$(mktemp -d)
OUT=$REL/dev/review-rel081-hunks.txt
: > "$OUT"
check() {
  f=$1; shift
  git -C "$REL" show "e0cce906:$f" 2>/dev/null | tr -d '\r' > "$T/base"
  tr -d '\r' < "$REL/$f" > "$T/rel"
  cp "$T/base" "$T/merged"
  for l in "$@"; do
    tr -d '\r' < "$R/frmtmb-wt-$l/$f" > "$T/$l"
    diff "$T/base" "$T/$l" | grep '^> ' | cut -c3- | grep -v '^[[:space:]]*$' > "$T/add"
    diff "$T/base" "$T/$l" | grep '^< ' | cut -c3- | grep -v '^[[:space:]]*$' > "$T/del"
    na=$(wc -l < "$T/add"); nd=$(wc -l < "$T/del")
    miss=0; kept=0
    while IFS= read -r line; do
      grep -qxF -- "$line" "$T/rel" || { miss=$((miss+1)); echo "  MISSING+ $l: ${line:0:150}" >> "$OUT"; }
    done < "$T/add"
    while IFS= read -r line; do
      if grep -qxF -- "$line" "$T/rel"; then
        nb=$(grep -cxF -- "$line" "$T/base"); nr=$(grep -cxF -- "$line" "$T/rel")
        [ "$nr" -lt "$nb" ] || { kept=$((kept+1)); echo "  KEPT- $l: ${line:0:150}" >> "$OUT"; }
      fi
    done < "$T/del"
    echo "$f $l: added $na missing $miss; removed $nd still-present $kept" >> "$OUT"
    git merge-file -q -p "$T/merged" "$T/base" "$T/$l" > "$T/m2"; rc=$?
    [ $rc -eq 0 ] || echo "  CONFLICT merging $l into $f ($rc)" >> "$OUT"
    cp "$T/m2" "$T/merged"
  done
  if cmp -s "$T/merged" "$T/rel"; then
    echo "$f: fresh merge EQUALS release" >> "$OUT"
  else
    echo "$f: fresh merge DIFFERS from release ($(diff "$T/merged" "$T/rel" | grep -c '^[<>]') lines)" >> "$OUT"
    diff "$T/merged" "$T/rel" | head -20 | cut -c1-160 | sed 's/^/    /' >> "$OUT"
  fi
}
check extensions/frmtmb.sample/R/loo.R hollow support
check dev/equiv-flow.R translator compat
check dev/equiv-specs.R translator compat
check dev/equiv-summary.R translator compat
check dev/equiv-findings.md translator compat
check dev/brmsport-verdicts-manual-fit.tsv code8 hollow
check dev/test-backlog.md compat code8 hollow
check NEWS.md compat code8 hollow support
check extensions/frmtmb.sample/NEWS.md hollow support
rm -rf "$T"
echo HUNKS DONE >> "$OUT"
