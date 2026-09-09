# Reviewer's private library and package attach. Separate from the
# lane's own dev/lincmt/lincmt-src.R so that the two never contend.
.libPaths(c("C:/Users/adf44/source/r/rev-lincmt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(RTMB)
  library(frmtmb)
  library(frmtmb.ode)
})
cat("frmtmb.ode from:", dirname(system.file(package = "frmtmb.ode")),
    "\n")
cat("frmtmb.ode version:",
    as.character(utils::packageVersion("frmtmb.ode")), "\n")
cat("RTMBode version:",
    as.character(utils::packageVersion("RTMBode")), "\n")
