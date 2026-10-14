#!/usr/bin/env bash
# Release reviewer of 0.71.0: an independent textual reconstruction of
# the seven-lane merge, written without reading dev/rel071-merge.sh's
# logic: every tracked file a lane changed against d2ef6312 (generated
# files left out), from its HEAD (253b3ef7) blob, git merge-file per lane
# in the order dblbar, famsimplex, famuni, postfix, optstart, illse,
# remeasure against the d2ef6312 blob, LF normalized; then each result
# against the release tree. No git write; work in dev/relrev071-out/recon.
set -u
REL=/c/Users/adf44/source/r/frmtmb-wt-release
W=$REL/dev/relrev071-out/recon
OUT=$REL/dev/relrev071-log/recon.txt
rm -rf "$W"; mkdir -p "$W"
excl='^(man/.*\.Rd|NAMESPACE|R/links-brms\.R|extensions/[^/]*/man/.*\.Rd|extensions/[^/]*/NAMESPACE|extensions/frmtmb.sample/tests/testthat/helper-brms-suite.R)$'
LANES="dblbar famsimplex famuni postfix optstart illse remeasure"
lf() { tr -d '\r'; }
: > "$W/files"
for l in $LANES; do
  git -C "/c/Users/adf44/source/r/frmtmb-wt-$l" diff --name-only d2ef6312 |
    grep -Ev "$excl" | sed "s|^|$l |" >> "$W/files"
done
cd "$REL"
: > "$OUT"
nf=0; nconf=0; ndiff=0
for f in $(awk '{print $2}' "$W/files" | sort -u); do
  nf=$((nf + 1))
  k=$(echo "$f" | tr '/' '_')
  git show "HEAD:$f" 2>/dev/null | lf > "$W/$k.cur" || : > "$W/$k.cur"
  git show "d2ef6312:$f" 2>/dev/null | lf > "$W/$k.base" || : > "$W/$k.base"
  who=""
  for l in $LANES; do
    grep -qx "$l $f" "$W/files" || continue
    who="$who $l"
    lf < "/c/Users/adf44/source/r/frmtmb-wt-$l/$f" > "$W/$k.$l"
    if ! git merge-file -q -p "$W/$k.cur" "$W/$k.base" "$W/$k.$l" \
         > "$W/$k.next" 2>/dev/null; then
      echo "CONFLICT $f at $l" >> "$OUT"; nconf=$((nconf + 1))
    fi
    mv "$W/$k.next" "$W/$k.cur"
  done
  if [ -f "$f" ]; then lf < "$f" > "$W/$k.rel"; else : > "$W/$k.rel"; fi
  if cmp -s "$W/$k.cur" "$W/$k.rel"; then
    echo "SAME $f ($who )" >> "$OUT"
  else
    n=$(diff "$W/$k.cur" "$W/$k.rel" | grep -c '^[<>]')
    echo "DIFF $f ($who ) $n lines" >> "$OUT"; ndiff=$((ndiff + 1))
    diff "$W/$k.cur" "$W/$k.rel" > "$W/$k.diff"
  fi
done
echo "files $nf conflicts $nconf differing $ndiff" >> "$OUT"
tail -n 1 "$OUT"
