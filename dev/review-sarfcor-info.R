.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
r <- tryCatch(frmtmb_covstruct_info("nope"), error = function(e) conditionMessage(e)); print(r)
r <- tryCatch(frmtmb_covstruct_info("car"), error = function(e) conditionMessage(e)); print(r)
suppressMessages(library(frmtmb.spatial))
str(frmtmb_covstruct_info("car")); str(frmtmb_covstruct_info("spde"))
r <- tryCatch(frmtmb_covstruct_info("sar"), error = function(e) conditionMessage(e)); print(r)
