# Release review 0.80.0: the test proposed for MAJOR 1, run on a library.
# It must fail on rellib-r19 and pass with the fix (rev-rel-lib/mut1).
#
#   Rscript dev/review-rel-mm22-test.R <lib list>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(strsplit(lib, ";")[[1]], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
source(file.path("C:/Users/adf44/source/r/frmtmb-wt-rel080/tests/testthat",
                 "helper-warnings.R"))
caught <- function(expr) {
  conds <- list()
  value <- withCallingHandlers(expr,
    warning = function(w) {
      conds[[length(conds) + 1L]] <<- w
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      conds[[length(conds) + 1L]] <<- m
      invokeRestart("muffleMessage")
    })
  list(value = value, conds = conds)
}
r <- testthat::test_that("a probe winner on a ridge is not made silent by a second descent", {
  # binomial Michaelis-Menten, review of lane nlstart, fresh row mm 22:
  # the winner of the sign probe runs to K -> Inf along the ridge to the
  # linear limit, 4.6 below the maximum. Refit from there, a second
  # descent in the body's Jacobian basis slid further along the ridge and
  # the standard-error warning that says so was lost
  set.seed(22001)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                  xq = round(runif(n, 0, 5), 3),
                  xb = round(runif(n, 0, 3000), 1),
                  g = factor(sample(1:15, n, TRUE)))
  d$ug <- rnorm(15, 0, 0.4)[d$g]
  d$nt <- 10L
  d$y <- rbinom(n, 10, plogis(3 * d$xp / (1 + d$xp) - 2))
  fm <- bf(y | trials(nt) ~ Vm * xp / (K + xp), Vm ~ 1, K ~ 1, nl = TRUE)
  got <- caught(frm(fm, family = binomial(), data = d))
  at <- caught(frm(fm, family = binomial(), data = d, se = FALSE,
                   start = list(beta = c("Vm_(Intercept)" = -0.1,
                                         "K_(Intercept)" = -6))))
  l_def <- as.numeric(logLik(got$value))
  l_max <- as.numeric(logLik(at$value))
  warned <- any(vapply(got$conds, inherits, NA, "warning"))
  # either the default fit reaches the maximum, or it says it may not
  cat("  l_def", l_def, " l_max", l_max, " warned", warned, "
")
  expect_true(l_max - l_def <= 1e-6 * abs(l_max) || warned)
})
cat("TEST RESULT", if (isTRUE(r)) "pass" else "fail", "\n")
