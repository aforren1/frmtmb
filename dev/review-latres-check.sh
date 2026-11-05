#!/usr/bin/env bash
# R CMD check --as-cran on core, tests left to the gated run, to read the NOTE.
source dev/latres-env.sh
cd "$LATRES_ROOT" || exit 1
out=dev/review-latres-check
mkdir -p "$out"
export PANDOC_DIR="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="$PANDOC_DIR:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_LIBS="$LATRES_LIB;$BASE_LIB;$USER_LIB"
cd "$out" || exit 1
"$RCMD" CMD build "$LATRES_ROOT" > build.log 2>&1 || { tail -20 build.log; exit 1; }
tb=$(ls -t frmtmb_*.tar.gz | head -1)
"$RCMD" CMD check --as-cran --no-tests "$tb" > check.log 2>&1
grep -E "^Status|NOTE$|WARNING$|ERROR$" check.log | head -20
