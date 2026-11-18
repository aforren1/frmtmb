# Lane nlband, fix round 1: run a test file on the lane's build with the
# two fixes of the round taken out, to see the new tests fail on the
# build the review read: nl_hole_located() always TRUE (no test that the
# data place the pole in the hole) and nl_cold_par() NULL (hole runs
# start the rest of the model at the fit only). The boundary-message
# filter in nl_fit_probed() stays in.
# Rscript dev/nlband-mockrun.R <package> <test file>
local({
  a <- commandArgs(trailingOnly = TRUE)
  Sys.setenv(FRMTMB_RELLIB = "C:/Users/adf44/source/r/wt-nlband-lib")
  .libPaths(c("C:/Users/adf44/source/r/wt-nlband-lib",
              "C:/Users/adf44/source/r/rellib-r22",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressMessages(library(frmtmb))
  ns <- asNamespace("frmtmb")
  for (nm in c("nl_hole_located", "nl_cold_par")) {
    f <- if (nm == "nl_hole_located") function(fit, ck) TRUE else
      function(frame, pn) NULL
    environment(f) <- ns
    unlockBinding(nm, ns)
    assign(nm, f, envir = ns)
    lockBinding(nm, ns)
  }
  cat("MOCK: nl_hole_located TRUE, nl_cold_par NULL\n")
})
source("dev/release/run-tests.R")
