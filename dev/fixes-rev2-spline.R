# Reviewer of lane fixes, re-check: frmtmb.spline's curve and derivative
# outputs under the s() basis switch. Writes one rds per lib; the
# companion block below compares them.
#   Rscript dev/fixes-rev2-spline.R <lib> <out.rds>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.spline)})
cat("LIB", find.package("frmtmb"), find.package("frmtmb.spline"), "\n")
set.seed(1)
dd <- data.frame(x = sort(runif(300)), g = factor(rep(1:15, 20)))
dd$y <- 2 * sin(pi * dd$x) + rnorm(15, 0, 0.5)[dd$g] + rnorm(300, 0, 0.4)
nd <- data.frame(x = seq(0.05, 0.95, length.out = 19))
fits <- list(
  s = frm(bf(y ~ s(x, k = 8)), data = dd),
  sre = frm(bf(y ~ s(x, k = 8) + (1 | g)), data = dd),
  sig = frm(bf(y ~ s(x, k = 8), sigma ~ s(x, k = 5)), data = dd))
out <- lapply(fits, function(f) {
  cv <- frm_curve(f, newdata = nd, nsim = 2000, seed = 7)
  dv <- frm_curve_deriv(f, var = "x", newdata = nd, nsim = 2000, seed = 7)
  ft <- frm_curve_feature(f, var = "x", type = "maximum",
                          newdata = data.frame(x = seq(0, 1, length.out = 200)))
  list(ll = as.numeric(logLik(f)), est = cv$.estimate, se = cv$.se,
       lsim = cv$.lower_sim, dest = dv[[grep("estimate", names(dv))[1]]],
       dse = dv[[grep("^[.]se$", names(dv))[1]]],
       feat = unlist(ft[1, vapply(ft, is.numeric, NA)]),
       fe = fixef(f)[, "Estimate"])
})
saveRDS(out, a[2])
