# Mutation check for item 3.2: with the epoch boundaries ignored (the
# spans found over the whole concatenated record, as if the epochs were
# one recording), test-epochs.R must fail. Runs the lane build with one
# internal function replaced in the namespace.
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
Sys.setenv(NOT_CRAN = "true")
suppressMessages({
  library(testthat); library(frmtmb); library(frmtmb.coupling)
})
setwd("C:/Users/adf44/source/r/frmtmb-wt-phase3a/extensions/frmtmb.coupling")
ns <- asNamespace("frmtmb.coupling")
orig <- get("cp_xspec_record", ns)
run <- function(label, fun) {
  unlockBinding("cp_xspec_record", ns)
  assign("cp_xspec_record", fun, envir = ns)
  lockBinding("cp_xspec_record", ns)
  res <- as.data.frame(testthat::test_file(
    "tests/testthat/test-epochs.R",
    reporter = testthat::SilentReporter$new(),
    env = testthat::test_env("frmtmb.coupling")))
  bad <- res$failed > 0 | res$error
  cat(sprintf("%-26s %d of %d blocks fail: %s\n", label, sum(bad),
              nrow(res), paste(substr(res$test[bad], 1, 45),
                               collapse = "; ")))
}
# the mutant: one epoch, the concatenation, so segments may cross
joined <- function(ex, ey, ...) {
  orig(list(unlist(ex, use.names = FALSE)),
       list(unlist(ey, use.names = FALSE)), ...)
}
run("unmutated", orig)
run("epoch boundaries ignored", joined)
