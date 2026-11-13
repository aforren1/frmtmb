# Lane boundfam, consolidation: the separated bernoulli of test-se-limit.R
# (seed 514), fitted() at the default fit and at a fit stopped short (the
# route OpenBLAS 0.3.26 takes at code 1): NaN count and warnings.
#   Rscript dev/boundfam-sepwarn.R <lib or "base">
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(514)
ds <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
ds$yb <- as.integer(ds$z > 0)
short <- frmtmb_control(optCtrl = list(eval.max = 40, iter.max = 30),
                        restarts = 0)
for (k in c("default", "short")) {
  f <- suppressWarnings(frm(yb ~ z + x, family = bernoulli(), data = ds,
                            control = if (k == "short") short else frmtmb_control()))
  w <- character()
  fe <- withCallingHandlers(fitted(f), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  lp <- withCallingHandlers(frm_linpred(f, se.fit = TRUE), warning = function(x) {
    w <<- c(w, paste("linpred:", conditionMessage(x))); invokeRestart("muffleWarning")
  })
  cat(sprintf("%s: code %d, se_explained %s, fitted NaN %d of %d, warnings %d\n",
              k, f$opt$convergence, format(f$cache$se_explained),
              sum(is.nan(fe[, "Est.Error"])), nrow(fe), length(w)))
  for (x in w) cat("   W:", substr(x, 1, 110), "\n")
}
