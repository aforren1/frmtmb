## The names the parser protects with a `.frm_` alias: cs_specials() is
## cs_names() less smooth, us_t and diag_t, plus rr and propto (R/parse.R).
##   Rscript dev/brmsnames-1024-probe0.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r12", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
print(frmtmb:::cs_names())
