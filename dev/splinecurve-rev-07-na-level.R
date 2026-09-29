# Reviewer, lane splinecurve: was an unseen bar-term level reachable
# through frm_curve() BEFORE the flag, by an NA grouping value?
#   Rscript splinecurve-rev-07-na-level.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spline)
})
cat("frmtmb.spline from", dirname(find.package("frmtmb.spline")), "\n")
set.seed(3)
dg <- data.frame(x = runif(300), g = factor(rep(1:15, each = 20)))
dg$y <- sin(2 * pi * dg$x) + rnorm(15, 0, 0.5)[dg$g] + rnorm(300, 0, 0.3)
fg <- frm(bf(y ~ s(x, k = 9) + (1 | g)), family = gaussian(), data = dg)
gna <- data.frame(x = seq(0, 1, length.out = 15),
                  g = factor(NA, levels = levels(dg$g)))
r <- tryCatch(frm_curve(fg, newdata = gna, re_formula = NULL, seed = 1,
                        nsim = 5000),
              error = function(e) e)
if (inherits(r, "error")) {
  cat("NA level, re_formula = NULL, default flag: REFUSED:",
      substr(conditionMessage(r), 1, 150), "\n")
} else {
  lb <- frm_lp_basis(fg, newdata = gna, re_formula = NULL)
  cat("NA level, re_formula = NULL, default flag: ANSWERED; extra_var",
      format(range(lb$extra_var), digits = 5), " crit_sim",
      format(r$.crit_sim[1], digits = 5), "\n")
}
