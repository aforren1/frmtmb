# Reviewer, final check: the hand-merged difference route of frmtmb.spline
# (the release's stacked grid plus the lane's lost rows) and
# frm_curve_feature(), on y ~ a + b with a, b ~ 1 + x (only a + b is
# determined). A difference of `a` between two x grids moves along a_x,
# which is lost; a difference of mu is determined.
#   Rscript dev/nanse-rev3-spline-diff.R merge|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.spline)})
cat("arm", arm, "spline", find.package("frmtmb.spline"),
    as.character(packageVersion("frmtmb.spline")), "\n")
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) {
    structure(conditionMessage(e), class = "err")
  }), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
set.seed(955)
dd <- data.frame(x = rnorm(60))
dd$y <- 3 + 0.5 * dd$x + rnorm(60, 0, 0.4)
fit <- suppressWarnings(frm(bf(y ~ a + b, a ~ 1 + x, b ~ 1 + x, nl = TRUE),
                            data = dd))
A <- data.frame(x = seq(-2, 2, length.out = 9))
B <- data.frame(x = seq(-1, 3, length.out = 9))
for (dp in list(NULL, "a")) {
  r <- cap(frm_curve(fit, newdata = A, contrast = B, dpar = dp,
                     simultaneous = FALSE))
  lab <- if (is.null(dp)) "mu" else dp
  if (inherits(r$v, "err")) {
    cat("difference", lab, "ERROR", substr(r$v, 1, 150), "\n")
  } else {
    cat(sprintf("difference %-2s: .se %s, finite %d of 9, warnings %d\n", lab,
                paste(signif(range(r$v$.se), 3), collapse = " to "),
                sum(is.finite(r$v$.se)), length(r$w)))
  }
}
# a feature (crossing) of a quadratic in x for a curve that crosses 0
set.seed(3)
d2 <- data.frame(x = runif(200, -3, 3))
d2$y <- d2$x^2 - 1 + rnorm(200, 0, 0.3)
f2 <- suppressWarnings(frm(bf(y ~ a + b, a ~ 1 + x + I(x^2), b ~ 1, nl = TRUE),
                           data = d2))
for (dp in list(NULL, "a")) {
  r <- cap(frm_curve_feature(f2, var = "x", type = "crossing", at = 0,
                             newdata = data.frame(x = seq(-3, 3, 0.1)),
                             dpar = dp))
  lab <- if (is.null(dp)) "mu" else dp
  if (inherits(r$v, "err")) {
    cat("feature", lab, "ERROR", substr(r$v, 1, 150), "\n")
  } else {
    cat(sprintf("feature %-2s: roots %s .se %s .value_se %s warnings %d\n", lab,
                paste(signif(r$v$.estimate, 3), collapse = " "),
                paste(signif(r$v$.se, 3), collapse = " "),
                paste(signif(r$v$.value_se, 3), collapse = " "), length(r$w)))
  }
}
