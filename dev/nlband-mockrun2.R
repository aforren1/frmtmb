# Lane nlband, fix round 2: run a test file on the lane's build with the
# round's two rules put back as round 1 had them, to see the new tests
# fail: nl_refit_held() with every message of the refit suppressed, and
# the located bound at half the hole's width (nl_hole_se() returns the
# standard error times span / half width, so the round-2 comparison with
# the span is round 1's with the half width). The ordered candidates and
# the refusal message of round 2 stay in.
# Rscript dev/nlband-mockrun2.R <package> <test file>
local({
  Sys.setenv(FRMTMB_RELLIB = "C:/Users/adf44/source/r/wt-nlband-lib")
  .libPaths(c("C:/Users/adf44/source/r/wt-nlband-lib",
              "C:/Users/adf44/source/r/rellib-r22",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressMessages(library(frmtmb))
  ns <- asNamespace("frmtmb")
  se_new <- ns$nl_hole_se
  f1 <- function(run, st) {
    held <- new.env(parent = emptyenv())
    held$w <- list()
    r <- tryCatch(
      withCallingHandlers(suppressMessages(run(st, FALSE)),
                          warning = function(w) {
                            held$w[[length(held$w) + 1L]] <- w
                            invokeRestart("muffleWarning")
                          }),
      error = function(e) NULL)
    list(fit = r, held = held$w)
  }
  f2 <- function(fit, ck) {
    se_new(fit, ck) * ck$span / ((ck$hi - ck$lo) / 2)
  }
  for (nm in c("nl_refit_held", "nl_hole_se")) {
    f <- if (nm == "nl_refit_held") f1 else f2
    # f2 keeps its own enclosure, which holds the real nl_hole_se()
    if (nm == "nl_refit_held") environment(f) <- ns
    unlockBinding(nm, ns)
    assign(nm, f, envir = ns)
    lockBinding(nm, ns)
  }
  cat("MOCK: round-1 nl_refit_held, located bound at half the width\n")
})
source("dev/release/run-tests.R")
