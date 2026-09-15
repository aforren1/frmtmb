export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export R_LIBS_USER="C:/Users/adf44/source/r/genrev2-lib;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=true
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-generics/dev/stan-cache
cd /c/Users/adf44/source/r/genrev2-check
rm -f frmtmb_*.tar.gz
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD build /c/Users/adf44/source/r/frmtmb-wt-generics
TB=$(ls -t frmtmb_*.tar.gz | head -1)
"/c/Program Files/R/R-4.6.1/bin/R.exe" CMD check --as-cran "$TB"
