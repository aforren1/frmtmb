# Which arm of test-extension-home-mixed.R's "family alone" block warns
# in predict(), and whether predict() of a fit whose family's package is
# absent refuses or answers.
#   Rscript dev/smallfix-predict-probe.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(testthat)})
env <- new.env()
sys.source("tests/testthat/helper-registry.R", envir = env)
f <- function(fam_pkg) {
  env$local_frmtmb_registries()
  frmtmb_register_covstruct("mxiid", lpdf = function(b, theta, blk) {
    sum(RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE))
  }, vcov = function(theta, blk) {
    matrix(exp(2 * theta[1]), 1, 1,
           dimnames = list(blk[["cnms"]], blk[["cnms"]]))
  }, npar = 1L, sd_idx = 1L, package = "frmtmbabsenttestpkgb")
  set.seed(1)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = stats::rnorm(100))
  d$y <- 1 + d$x + stats::rnorm(20)[d$g] + stats::rnorm(100, 0, 0.5)
  fit <- frm(bf(y ~ x + mxiid(1 | g)), family = gaussian(), data = d)
  attr(fit$spec$responses[[1L]]$family, "frmtmb_package") <- fam_pkg
  fit$frame$spec <- fit$spec
  fit$cache <- new.env()
  for (m in c("summary", "vcov", "predict")) {
    r <- withCallingHandlers(
      tryCatch({ get(m)(fit); "ok" }, error = function(e)
        paste("ERROR", class(e)[1], substr(conditionMessage(e), 1, 90))),
      warning = function(w) {
        cat("   WARNING in", m, ":", substr(conditionMessage(w), 1, 90), "\n")
        invokeRestart("muffleWarning")
      })
    cat(fam_pkg, m, ":", r, "\n")
  }
}
f("frmtmbabsenttestpkga")
f("frmtmb")
