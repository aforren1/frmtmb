#!/usr/bin/env bash
# R CMD build and check --as-cran of one package inside dev/gpgap-check/.
#   bash dev/gpgap-check.sh <package dir relative to the worktree>
set -u
root="$(pwd)"
pkgdir="$1"
name="$(basename "$root/$pkgdir")"
[ "$pkgdir" = "." ] && name=frmtmb
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
rt=/c/rtools45
pan="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
tex=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
export PATH="$pan:$tex:$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_LIBS="C:/Users/adf44/source/r/wt-gpgap-lib;C:/Users/adf44/source/r/rellib-r16;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
out="dev/gpgap-check/$name"
rm -rf "$out"; mkdir -p "$out"
cd "$out"
"$R" CMD build "$root/$pkgdir" > build.txt 2>&1
tb=$(ls ./*.tar.gz | head -1)
"$R" CMD check --as-cran "$tb" > check.txt 2>&1
tail -3 check.txt
