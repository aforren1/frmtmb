# Reviewer, lane splinecurve, re-check: the refusal's text when the
# nonzero extra_var is an exact gp() kriging variance, with no unseen
# level anywhere.
#   Rscript splinecurve-rev-10-gpmsg.R > splinecurve-rev-log/r2-10-gpmsg.txt
LIB <- "C:/Users/adf44/source/r/wt-splinecurve-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
set.seed(8)
dgp <- data.frame(x = sort(runif(60, 0, 5)))
dgp$y <- sin(dgp$x) + rnorm(60, 0, 0.3)
fgp <- frm(bf(y ~ gp(x)), family = gaussian(), data = dgp)
gx <- data.frame(x = dgp$x[-1] - diff(dgp$x) / 2)
e <- tryCatch(frm_curve_deriv(fgp, var = "x", newdata = gx,
                              allow_new_levels = TRUE, simultaneous = FALSE),
              error = function(e) e)
cat(strwrap(conditionMessage(e), 78), sep = "\n")
r <- tryCatch(frm_curve_deriv(fgp, var = "x", newdata = gx, re_formula = NA,
                              allow_new_levels = TRUE, simultaneous = FALSE),
              error = function(e) e)
cat("\nfollowing the advice, re_formula = NA:",
    if (inherits(r, "error")) "still REFUSED" else "answered", "\n")
