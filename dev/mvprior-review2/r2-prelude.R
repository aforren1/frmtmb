# Reviewer 2 prelude. ARM=lane reads the lane library first; ARM=base
# reads rellib-r3 first. User library last.
r2_arm <- Sys.getenv("R2_ARM", "lane")
R2_ROOT <- "C:/Users/adf44/source/r/frmtmb-wt-mvprior"
.libPaths(c(if (identical(r2_arm, "lane")) "C:/Users/adf44/source/r/mvprior-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = file.path(R2_ROOT, "dev/stan-cache"),
           NOT_CRAN = "true")
suppressPackageStartupMessages(library(frmtmb))
cat("ARM", r2_arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")

# evaluate an expression, returning "OK" or the error message
try_msg <- function(expr) {
  r <- tryCatch({force(expr); "OK"}, error = function(e) {
    paste("ERR:", gsub("\\s+", " ", conditionMessage(e)))
  })
  r
}
