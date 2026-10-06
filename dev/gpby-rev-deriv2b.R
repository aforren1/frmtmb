# Reviewer: frm_curve_deriv(order = 2) on an exact gp() past its data,
# both arms, and the lane's frm_curve_deriv(order = 2) with a larger eps.
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spline"), "\n")
set.seed(5)
x <- sort(runif(50, 0, 4))
dg <- data.frame(x = x, y = sin(1.3 * x) + rnorm(50, 0, 0.15))
fg <- frm(bf(y ~ gp(x)), data = dg)
xs <- c(3.0, 4.2, 4.8, 5.5)
g2 <- frm_curve_deriv(fg, var = "x", order = 2, newdata = data.frame(x = xs),
                      simultaneous = FALSE)
cat("order 2 se:", sprintf("%.6g", g2$.se), "| eps", attr(g2, "eps"), "\n")
g2b <- frm_curve_deriv(fg, var = "x", order = 2, newdata = data.frame(x = xs),
                       simultaneous = FALSE, eps = 6e-3)
cat("order 2 se, eps 6e-3:", sprintf("%.6g", g2b$.se), "\n")
g2s <- tryCatch(frm_curve_deriv(fg, var = "x", order = 2,
                                newdata = data.frame(x = seq(3.5, 5.5, length.out = 21)),
                                simultaneous = TRUE, nsim = 2000, seed = 1),
                error = function(e) e)
if (inherits(g2s, "error")) cat("simultaneous order 2:", conditionMessage(g2s), "\n") else
  cat("simultaneous order 2 crit", g2s$.crit_sim[1], "se", sprintf("%.4g", g2s$.se[c(1, 6, 11, 16, 21)]), "\n")
g1 <- frm_curve_deriv(fg, var = "x", order = 1, newdata = data.frame(x = xs),
                      simultaneous = FALSE)
cat("order 1 se:", sprintf("%.6g", g1$.se), "\n")
