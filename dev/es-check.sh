#!/bin/sh
ES_LIB="C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/es-lib"
OUT="C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/es-check"
WT="C:/Users/adf44/source/r/frmtmb-wt-esicar"
mkdir -p "$OUT"
cd "$OUT" || exit 1
export R_LIBS="$ES_LIB;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export RSTUDIO_PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="$RSTUDIO_PANDOC:$PATH"
export _R_CHECK_CRAN_INCOMING_=false
export NOT_CRAN=false
echo "pandoc: $(command -v pandoc)"
"C:/Program Files/R/R-4.6.1/bin/R.exe" CMD build --no-build-vignettes "$WT" || exit 1
TB=$(ls -t frmtmb_*.tar.gz | head -1)
echo "tarball: $TB"
"C:/Program Files/R/R-4.6.1/bin/R.exe" CMD check --as-cran --no-manual \
  --library="$ES_LIB" -o "$OUT" "$TB"
echo "CHECK EXIT $?"
