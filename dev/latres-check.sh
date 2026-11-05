#!/usr/bin/env bash
# R CMD check --as-cran once per changed package, inside dev/latres-check/.
#   bash dev/latres-check.sh <package dir relative to the worktree> <name>
source "$(dirname "$0")/latres-env.sh"
cd "$LATRES_ROOT" || exit 1
pkgdir="$1"; name="$2"
out="dev/latres-check/$name"
mkdir -p "$out"
export PANDOC_DIR="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="$PANDOC_DIR:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_LIBS="$LATRES_LIB;$BASE_LIB;$USER_LIB"
abs="$LATRES_ROOT/$pkgdir"
[ "$pkgdir" = "." ] && abs="$LATRES_ROOT"
cd "$out" || exit 1
"$RCMD" CMD build "$abs" > build.log 2>&1 || { tail -20 build.log; exit 1; }
tb=$(ls -t "$name"_*.tar.gz | head -1)
"$RCMD" CMD check --as-cran "$tb" > check.log 2>&1
grep -E "^Status|^\* DONE|^ERROR|WARNING$" check.log | head -20
tail -3 check.log | cut -c1-160
