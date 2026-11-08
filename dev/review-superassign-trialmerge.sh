#!/usr/bin/env bash
# Review of lane superassign: pairwise trial merges, superassign against
# each other lane, of every R/ or extensions/ file both change.
# git merge-file on LF-normalized copies (base = main checkout,
# 3cd2a0c8). Output tree: C:/Users/adf44/source/r/rev-superassign-src/merge.
# Prints conflicts per file and parses every merged file for `<<-`.
set -u
R=C:/Users/adf44/source/r
OUT=$R/rev-superassign-src/merge
rm -rf "$OUT"; mkdir -p "$OUT"
SA=$R/frmtmb-wt-superassign
export GIT_OPTIONAL_LOCKS=0
files_of() { git -C "$1" diff --name-only -- R extensions 2>/dev/null; }
for L in fitedge nlstart readers loglik nestquad arstep bigx osa latperf mfxse ci; do
  W=$R/frmtmb-wt-$L
  for F in $(comm -12 <(files_of "$SA" | sort) <(files_of "$W" | sort)); do
    case "$F" in *.R) ;; *) continue ;; esac
    d=$OUT/$L/$(dirname "$F"); mkdir -p "$d"
    b=$d/$(basename "$F").base; o=$d/$(basename "$F").sa; t=$d/$(basename "$F").$L
    tr -d '\r' < "$R/frmtmb/$F" > "$b"
    tr -d '\r' < "$SA/$F" > "$o"
    tr -d '\r' < "$W/$F" > "$t"
    m=$d/$(basename "$F")
    cp "$o" "$m"
    git merge-file -L superassign -L base -L "$L" "$m" "$b" "$t"
    n=$?
    printf 'MERGE %-9s %-34s conflicts %d\n' "$L" "$F" "$n"
  done
done
