#!/usr/bin/env bash
# Lane seflags, punch round 1: a trial merge of this lane with lane
# optrap (its par_units for identity-link means), built in
# dev/seflags-merge-out/src and installed into the scratch library
# dev/seflags-merge-out/lib (both gitignored). optrap's tree is read,
# never written: its R files that this lane does not touch, and its
# changed test files, are copied over a copy of this lane's package.
#   bash dev/seflags-merge.sh
set -eu
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
OP=/c/Users/adf44/source/r/frmtmb-wt-optrap
O=$W/dev/seflags-merge-out
rm -rf "$O/src"
mkdir -p "$O/src" "$O/lib"
cd "$W"
for x in DESCRIPTION NAMESPACE R man tests inst data vignettes .Rbuildignore; do
  [ -e "$x" ] && cp -r "$x" "$O/src/"
done
# the files optrap changed (this lane changes none of them)
mine=$(git -C "$W" diff --name-only)
for f in $(git -C "$OP" diff --name-only); do
  case "$f" in R/*|tests/*)
    if echo "$mine" | grep -qx "$f"; then echo "OVERLAP $f"; exit 1; fi
    cp "$OP/$f" "$O/src/$f"; echo "from optrap: $f" ;;
  esac
done
for f in $(git -C "$OP" ls-files --others --exclude-standard R tests); do
  cp "$OP/$f" "$O/src/$f"; echo "from optrap (new): $f"
done
echo "optrap HEAD $(git -C "$OP" rev-parse --short HEAD), copied $(date)"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL \
  --library="C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-merge-out/lib" \
  --no-multiarch "C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-merge-out/src"
