# Lane wt-defects: the preamble every probe sources. The library arm is
# chosen by DEFECTS_ARM: "before" reads rellib-r3 (base 0.65.0), anything
# else puts the lane's private library first.
arm <- Sys.getenv("DEFECTS_ARM", "after")
libs <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (!identical(arm, "before")) {
  libs <- c("C:/Users/adf44/source/r/wt-defects-lib", libs)
}
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = "frmtmb")
suppressMessages({
  library(testthat)
  library(frmtmb)
})
cat("ARM", arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
source("tests/testthat/helper-brms-suite.R")
# print a condition or value compactly
show <- function(expr) {
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    invisible(e)
  })
  r
}
