#!/usr/bin/env bash
# Fix round 1: three-way trial of ncgp against each other lane over the base, per file both touch.
# Writes only under rev-ncgp-src/merge2. No git state is changed; `git status` is read only.
export PATH="/c/rtools45/usr/bin:$PATH"
cd /c/Users/adf44/source/r
S=rev-ncgp-src/merge2; rm -rf $S; mkdir -p $S
base=frmtmb-wt-r16base; a=frmtmb-wt-ncgp; T=rev-ncgp-src/t2
files() { (cd "$1"; git status --short | grep -v "^?? dev/\|^ M dev/\|^?? tests\|^?? extensions/.*/tests" | cut -c4-); }
files $a | sort > $T.a.list
for lane in gpgap refuse latres nlfix predread impsamp; do
  b=frmtmb-wt-$lane
  files $b | sort > $T.b.list
  both=$(comm -12 $T.a.list $T.b.list)
  tot=0
  for f in $both; do
    [ -f "$base/$f" ] || { echo "$lane $f new on both"; continue; }
    mkdir -p "$S/$lane/$(dirname $f)"
    tr -d '\r' < $a/$f > $T.a; tr -d '\r' < $base/$f > $T.o; tr -d '\r' < $b/$f > $T.b
    diff3 -m -L ncgp -L base -L $lane $T.a $T.o $T.b > "$S/$lane/$f"
    n=$(grep -c '^<<<<<<< ncgp' "$S/$lane/$f")
    tot=$((tot + n))
    echo "$lane $f conflicts=$n"
  done
  echo "TOTAL $lane $tot"
done
rm -f $T.a $T.o $T.b $T.a.list $T.b.list
