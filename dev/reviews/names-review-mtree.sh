#!/usr/bin/env bash
# Review of lane names: build names + sampfront merged (core and
# frmtmb.sample) in the scratchpad and install it into a scratch library.
# Run dev/reviews/names-review-merge.sh first (it writes the merged files).
set -u
export GIT_OPTIONAL_LOCKS=0
SP=/c/Users/adf44/AppData/Local/Temp/2/claude/ee0bfefdf562474c/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad
N=/c/Users/adf44/source/r/frmtmb-wt-names
S=/c/Users/adf44/source/r/frmtmb-wt-sampfront
M=$N/dev/reviews/names-review-merge
T=$SP/mtree
L=$SP/mlib
rm -rf "$T" "$L"; mkdir -p "$T/frmtmb" "$L"
for d in DESCRIPTION NAMESPACE R man inst data tests vignettes .Rbuildignore; do
  cp -r "$N/$d" "$T/frmtmb/"
done
cp -r "$N/extensions/frmtmb.sample" "$T/frmtmb.sample"
rm -rf "$T/frmtmb.sample/dev"
cd "$S" || exit 1
for f in $( { git diff --name-only e8abed0a; git ls-files --others --exclude-standard; } | grep -v '^dev/' | grep -v 'NEWS.md$' | sort -u); do
  case "$f" in
    extensions/frmtmb.sample/*) dst="$T/frmtmb.sample/${f#extensions/frmtmb.sample/}";;
    extensions/*) continue;;
    *) dst="$T/frmtmb/$f";;
  esac
  k=$(echo "$f" | tr '/' '_')
  if [ -f "$M/sampfront.$k.merged" ]; then
    cp "$M/sampfront.$k.merged" "$dst"; echo "merged $f"
  else
    mkdir -p "$(dirname "$dst")"; cp "$f" "$dst"; echo "theirs $f"
  fi
done
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R=/c/Program\ Files/R/R-4.6.1/bin/R.exe
export R_LIBS="C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
cd "$T"
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library="$L" --no-multiarch frmtmb > "$T/install-core.log" 2>&1; echo "core install $?"
export R_LIBS="$L;C:/Users/adf44/source/r/rellib-r20;C:/Users/adf44/AppData/Local/R/win-library/4.6"
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library="$L" --no-multiarch frmtmb.sample > "$T/install-sample.log" 2>&1; echo "sample install $?"
