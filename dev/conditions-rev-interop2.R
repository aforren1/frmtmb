# Reviewer, lane wt-conditions (round 2): insight's
# .get_predicted_ci_modelmatrix() tests inherits(mm, "simpleError") on
# what get_modelmatrix() raised. Reach it with a frmtmb refusal from
# model.matrix(): a multivariate fit.
#   Rscript dev/conditions-rev-interop2.R base|lane
arm <- commandArgs(TRUE)[1L]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/conditions-lib"
            else "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20260917)
d <- data.frame(y = rnorm(60), z = rnorm(60), x = rnorm(60))
mv <- frm(mvbf(bf(y ~ x), bf(z ~ x)), data = d)
r <- function(label, expr) {
  cat("\n##", label, "\n")
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      cat("ERROR [", paste(class(e), collapse = "/"), "]",
          conditionMessage(e), "\n"); NULL }),
    warning = function(w) { cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning") })
  if (!is.null(v)) cat(utils::capture.output(str(v, max.level = 1)),
                       sep = "\n")
}
r("model.matrix(mv)", model.matrix(mv))
r("insight::get_modelmatrix(mv, data = d)",
  insight::get_modelmatrix(mv, data = d))
r(".get_predicted_ci_modelmatrix(mv, data = d)",
  insight:::.get_predicted_ci_modelmatrix(mv, data = d))
r("insight::get_predicted(mv, data = d, ci = 0.95)",
  insight::get_predicted(mv, data = d, ci = 0.95))
p <- rep(0, nrow(d))
r("insight::get_predicted_ci(mv, predictions = p, data = d)",
  insight::get_predicted_ci(mv, predictions = p, data = d))
g <- frm(y ~ x, data = d, family = gaussian())
r("insight::get_predicted_ci(gaussian, data = d, dpar = 'nosuch')",
  insight::get_predicted_ci(g, predictions = p, data = d, dpar = "nosuch"))
d$yp <- 2 * exp(-0.5 * abs(d$x)) + rnorm(60, sd = 0.1)
nl <- frm(bf(yp ~ a * exp(-b * abs(x)), a ~ 1, b ~ 1, nl = TRUE), data = d,
          start = list(beta = c(a_Intercept = 2, b_Intercept = 0.5)))
r("model.matrix(nl)", model.matrix(nl))
r("insight::get_predicted_ci(nl, predictions = p, data = d)",
  insight::get_predicted_ci(nl, predictions = p, data = d))
r("insight::get_predicted(nl, ci = 0.95)",
  insight::get_predicted(nl, ci = 0.95))
