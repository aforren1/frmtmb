## The margin the shipped assertion has at the seeds it pins.
##
## The test in tests/testthat/test-coherence.R compares
## mean(ratio) in the treatment arm against max(ratio) in the null arm
## over seeds 2610:2613. This prints both arms so the margin is a
## recorded number rather than something a reader has to re-run.
##
## Run: Rscript dev/coh-shipseeds.R
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
setwd("extensions/frmtmb.coupling/tests/testthat")
## The test file's own generator and fitters, so that this measures the
## shipped code rather than a copy of it. test_that() is needed because
## the file calls it at the top level.
test_that <- function(...) invisible(NULL)
source("test-coherence.R", local = TRUE)

ratio <- function(sd_ic, seeds) {
  vapply(seeds, function(s) {
    d <- cp_cells(s, n_sub = 16L, n_rep = 6L, sd_idcond = sd_ic)
    cp_width("cond + (1 | id) + (1 | id:cond)", d) /
      cp_width("cond + (1 | id)", d)
  }, numeric(1))
}
seeds <- 2610:2613
tr <- ratio(0.5, seeds)
nu <- ratio(0, seeds)
cat("treat:", paste(signif(tr, 5), collapse = " "), "\n")
cat("null :", paste(signif(nu, 5), collapse = " "), "\n")
cat("mean(treat) =", mean(tr), " max(null) =", max(nu),
    " margin =", mean(tr) / max(nu), "\n")
cat("min(treat) =", min(tr), "\n")
