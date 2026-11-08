# Release review 0.80.0: the test proposed for MAJOR 2, run on a library.
# It must fail on rellib-r19 and pass with the fix (rev-rel-lib/mut2).
#
#   Rscript dev/review-rel-noint-test.R <lib list>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(strsplit(lib, ";")[[1]], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
r <- testthat::test_that("a product of shifted columns beside `0 + f` is judged centered", {
  # `0 + f` spans the constant with no intercept column, and the rank
  # test read such a design uncentered: it dropped w1:w2 at 1e5 (capped)
  # and at 1e6 with the message "rank deficient", 84 to 101
  # log-likelihood units short (release review of 0.80.0)
  set.seed(7)
  n <- 300
  d <- data.frame(z1 = stats::rnorm(n), z2 = stats::rnorm(n),
                  f = factor(sample(letters[1:3], n, TRUE)))
  d$y <- c(a = 0, b = 1, c = -1)[as.character(d$f)] + 0.5 * d$z1 -
    0.4 * d$z2 + 0.6 * d$z1 * d$z2 + stats::rnorm(n, 0, 0.5)
  fz <- frm(y ~ 0 + f + z1 * z2, data = d)
  for (m in c(1e5, 1e6)) {
    d$w1 <- d$z1 + m
    d$w2 <- d$z2 + m
    cls <- character(0)
    fw <- NULL
    expect_no_message(fw <- withCallingHandlers(
      frm(y ~ 0 + f + w1 * w2, data = d),
      warning = function(x) {
        cls <<- c(cls, if (inherits(x, "frmtmb_design_capped")) "capped"
                  else class(x)[1L])
        invokeRestart("muffleWarning")
      }))
    cat("  m", m, "beta", length(fw$estimates[["beta"]]), "warnings",
        paste(cls, collapse = ","), "rel", as.numeric(logLik(fw)) /
          as.numeric(logLik(fz)) - 1, "\n")
    if (m == 1e5) {
      # the centered design is well conditioned here: kept and fitted
      expect_length(fw$estimates[["beta"]], 6L)
      expect_lt(abs(as.numeric(logLik(fw)) / as.numeric(logLik(fz)) - 1),
                1e-6)
    } else if (length(fw$estimates[["beta"]]) == 6L) {
      expect_lt(abs(as.numeric(logLik(fw)) / as.numeric(logLik(fz)) - 1),
                1e-6)
    } else {
      expect_true("capped" %in% cls)
    }
  }
})
cat("TEST RESULT", if (isTRUE(r)) "pass" else "fail", "\n")
