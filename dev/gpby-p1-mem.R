# Punch round 1, m5: memory of extra_cov and of frm_curve() on a
# 2000-row grid, after extra_cov became sparse (or absent) where nothing
# contributes. Peak is gc()'s max Vcells since a reset; the baseline is
# what was in use at the reset, so the call's own peak is the difference.
# Usage: Rscript dev/gpby-p1-mem.R <lane|base>
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.spline)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spline"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
src <- parse(file.path(wt, "tests/testthat/test-gp-by.R"))
for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e)
d <- gpby_data()
fit <- frm(bf(y ~ gp(x)), data = d)
fs <- frm(bf(y ~ s(x, k = 10)), data = d)
big <- data.frame(x = seq(-1, 9, length.out = 2000))
bigs <- data.frame(x = seq(0.1, 5.9, length.out = 2000))
peak <- function(lab, expr) {
  g0 <- gc(reset = TRUE)
  t <- system.time(r <- eval.parent(substitute(expr)))
  g1 <- gc()
  cat(sprintf("MEM %s %-44s peak %6.1f MB, baseline %5.1f MB, own %6.1f MB, %.2f s, result %.1f MB\n",
              arm, lab, g1[2, 6], g0[2, 2], g1[2, 6] - g0[2, 2],
              t[["elapsed"]], as.numeric(object.size(r)) / 2^20))
  invisible(NULL)
}
has_ec <- "extra_cov" %in% names(formals(frm_lp_basis))
if (has_ec) {
  peak("gp(x) frm_lp_basis(extra_cov = TRUE)",
       frm_lp_basis(fit, newdata = big, extra_cov = TRUE))
  peak("s(x) frm_lp_basis(extra_cov = TRUE)",
       frm_lp_basis(fs, newdata = bigs, extra_cov = TRUE))
  E <- frm_lp_basis(fit, newdata = big, extra_cov = TRUE)$extra_cov
  cat("gp extra_cov class:", class(E)[1], "\n")
  E <- frm_lp_basis(fs, newdata = bigs, extra_cov = TRUE)$extra_cov
  cat("s(x) extra_cov class:", class(E)[1], "nnz", length(E@x), "\n")
}
peak("gp(x) frm_lp_basis()", frm_lp_basis(fit, newdata = big))
peak("gp(x) frm_curve(nsim = 1000)",
     frm_curve(fit, newdata = big, nsim = 1000, seed = 1))
peak("s(x) frm_curve(nsim = 1000)",
     frm_curve(fs, newdata = bigs, nsim = 1000, seed = 1))
peak("s(x) frm_curve_deriv(nsim = 1000)",
     frm_curve_deriv(fs, var = "x", newdata = bigs, nsim = 1000, seed = 1))
cat("DONE\n")
