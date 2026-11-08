#!/usr/bin/env bash
# Review of lane superassign: for every file superassign shares with two
# or more lanes, merge the other lanes first (in the order listed), then
# superassign last, and the reverse (superassign first). Counts the
# conflicts each step adds. LF-normalized copies, base = main checkout.
set -u
R=C:/Users/adf44/source/r
OUT=$R/rev-superassign-src/merge2
rm -rf "$OUT"; mkdir -p "$OUT"
export GIT_OPTIONAL_LOCKS=0
lanes="fitedge nlstart readers loglik nestquad arstep bigx osa latperf mfxse"
files_of() { git -C "$1" diff --name-only -- R extensions 2>/dev/null; }
for F in $(files_of $R/frmtmb-wt-superassign | grep '\.R$'); do
  users=""
  for L in $lanes; do
    files_of $R/frmtmb-wt-$L | grep -qx "$F" && users="$users $L"
  done
  [ -z "$users" ] && continue
  b=$OUT/$(basename "$F").base
  tr -d '\r' < "$R/frmtmb/$F" > "$b"
  for order in last first; do
    m=$OUT/$(basename "$F").$order
    cp "$b" "$m"
    seq="$users"
    [ $order = last ] && seq="$users superassign" || seq="superassign $users"
    msg=""
    for L in $seq; do
      t=$OUT/$(basename "$F").$L
      tr -d '\r' < "$R/frmtmb-wt-$L/$F" > "$t"
      git merge-file -q "$m" "$b" "$t"
      n=$?
      msg="$msg $L:$n"
    done
    printf 'SEQ %-22s superassign-%-5s%s\n' "$F" "$order" "$msg"
  done
done
