# With core's gp() key back on brms's 15-digit rule, does frmtmb.spline's
# exact row key (sp_row_key()) disagree with core in a way that matters?
# Core merges two unseen gp() rows equal to 15 digits into one kriged
# position, but every other part of a row's design (a smooth, a linear
# term) is still computed from that row's own x. So at u and
# u * (1 + 2^-52) core's linear predictor can differ, and a spline key
# that merged the two rows would report row 1's value for row 2. This
# prints core's two values and the spline's, both arms of the key.
# Usage: Rscript dev/ciharden-splinekey.R <lib or base>
a <- commandArgs(TRUE)
.libPaths(c(if (!identical(a[1], "base")) a[1],
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(17)
xg <- round(runif(120, 0, 10), 1)
dg <- data.frame(y = sin(xg) + 0.2 * xg + rnorm(120, 0, 0.3), x = xg)
fit <- suppressWarnings(frm(bf(y ~ x + gp(x)), data = dg))
for (u in c(10 / 3, 7 / 3, 2.2 / 3, 5.123456789)) {
  nd <- data.frame(x = c(u, u * (1 + 2^-52)))
  p <- frm_linpred(fit, newdata = nd, se.fit = TRUE)
  lb <- frm_lp_basis(fit, newdata = nd)
  A <- unname(as.matrix(lb$A))
  cat(sprintf("u = %.10g: x identical %s; core eta %a %a (equal %s);",
              u, identical(nd$x[1], nd$x[2]), p$fit[1], p$fit[2],
              identical(unname(p$fit[1]), unname(p$fit[2]))),
      "gp columns equal:", identical(A[1, -(1:2)], A[2, -(1:2)]),
      "; x column equal:", identical(A[1, 2], A[2, 2]),
      "; extra_var equal:", identical(lb$extra_var[1], lb$extra_var[2]),
      "\n")
  k_exact <- frmtmb.spline:::sp_row_key(nd)
  k_paste <- do.call(paste, c(nd, sep = "\r"))
  cat("   spline exact key merges:", k_exact[1] == k_exact[2],
      "; a 15-digit key would merge:", k_paste[1] == k_paste[2], "\n")
}
