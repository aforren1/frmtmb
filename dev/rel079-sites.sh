#!/usr/bin/env bash
# The site checks of the 0.79.0 consolidation, as CI runs them:
#   docs   the nine pkgdown sites, dev/release/build-docs.R, into the
#          session scratchpad (not kept), R_LIBS rellib-r17 then the user
#          library, pandoc on PATH
#   saved  .github/scripts/saved-fit-check.R on a SCRATCH COPY of
#          rellib-r17 in the session scratchpad (the script moves
#          frmtmb.survival and frmtmb.spatial out of the first library
#          for the test and back); the copy is deleted after
#   wf     every workflow parses (dev/rel079-workflows.R)
# Adapted from the 0.78.0 record's "Site checks".
#   bash dev/rel079-sites.sh <docs|saved|wf>
set -u
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel079
SP=C:/Users/adf44/AppData/Local/Temp/2/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
cd "$ROOT" || exit 1
case "$1" in
  docs)
    s=$(date +%s)
    R_LIBS="C:/Users/adf44/source/r/rellib-r17;$USR" "$RS" dev/release/build-docs.R \
      --dest="$SP/pkgdown-site" --min-sites=9
    echo "exit $? after $(( $(date +%s) - s )) s"
    rm -rf "$SP/pkgdown-site";;
  saved)
    C="$SP/rellib-r17-copy"
    rm -rf "$C"; cp -r C:/Users/adf44/source/r/rellib-r17 "$C"
    ls "$C" | tr '\n' ' '; echo
    R_LIBS="$C;$USR" "$RS" .github/scripts/saved-fit-check.R
    echo "exit $?"
    echo "after: $(ls "$C" | tr '\n' ' ')"
    ls -d "$SP"/rellib-r17-copy* 2>/dev/null
    rm -rf "$C";;
  wf) "$RS" dev/rel079-workflows.R;;
  *) echo "unknown $1"; exit 1;;
esac
