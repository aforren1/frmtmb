# Punch round 1, m5: how much of frm_curve()'s gc() peak on a 2000-row
# gp(x) grid is the extra covariance, and how much the run-to-run spread
# of the peak is. Usage: Rscript dev/gpby-p1-mem2.R <lane|base>
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
big <- data.frame(x = seq(-1, 9, length.out = 2000))
for (k in 1:3) {
  for (sim in c(0, 1000)) {
    invisible(gc())
    g0 <- gc(reset = TRUE)
    r <- if (sim) frm_curve(fit, newdata = big, nsim = sim, seed = 1) else
      frm_curve(fit, newdata = big, simultaneous = FALSE)
    g1 <- gc()
    cat(sprintf("MEM %s run %d gp(x) frm_curve(nsim = %4d) own peak %6.1f MB\n",
                arm, k, sim, g1[2, 6] - g0[2, 2]))
    rm(r)
  }
}
