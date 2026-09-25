# Reviewer prelude for lane wt-simnewdata. RV_LIB=base drops the lane
# library so the same probe measures the 0.62.0 base build (rellib-r3).
lane_lib <- "C:/Users/adf44/source/r/simnewdata-lib"
rv_base <- identical(Sys.getenv("RV_LIB"), "base")
.libPaths(c(if (!rv_base) lane_lib,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-simnewdata/dev/stan-cache",
           NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", find.package("frmtmb"), "\n")
try_msg <- function(expr) {
  tryCatch({ force(expr); "ANSWERS" },
           error = function(e) paste("ERROR:", conditionMessage(e)))
}
