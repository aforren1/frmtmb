#!/usr/bin/env bash
# Lane wt-arcovsample: is the PDF-manual WARNING mine or the machine's?
# R CMD Rd2pdf with the index, on the BASE sources (main checkout,
# read-only) and on the LANE sources, same TinyTeX and same PATH. If the
# base sources warn too, the WARNING is environmental.
#
#   bash dev/arcovsample-rd2pdf.sh
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-arcovsample
RBIN="/c/Program Files/R/R-4.6.1/bin"
export R_LIBS="C:/Users/adf44/source/r/wt-arcovsample-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
O="$W/dev/arcovsample-log/rd2pdf"
rm -rf "$O"; mkdir -p "$O"
for arm in base lane; do
  case $arm in
    base) SRC="C:/Users/adf44/source/r/frmtmb/extensions/frmtmb.sample" ;;
    lane) SRC="$W/extensions/frmtmb.sample" ;;
  esac
  cd "$O"
  "$RBIN/R.exe" CMD Rd2pdf --no-preview --force \
    --output="$O/$arm.pdf" "$SRC" > "$O/$arm.log" 2>&1
  echo "$arm exit=$? pdf=$(ls -la "$O/$arm.pdf" 2>/dev/null | awk '{print $5}')"
  echo "  dest warnings: $(grep -c 'pdfTeX warning (dest)' "$O/$arm.log")"
  echo "  bang lines   : $(grep -c '^! ' "$O/$arm.log")"
done
