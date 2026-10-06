# Reviewer, final check: the release's predict() (predict-brms.R) draws
# the outer parameters from N(theta_hat, vcov(full = TRUE)). On a fit
# with a lost SE that matrix has NaN rows, so predict() falls back to
# the estimates and warns. Fixture: test-portability's singular
# (1 + x | g) block, where the fixed effects are determined. Compares
# Est.Error and the 95% interval width with the release build.
#   Rscript dev/nanse-rev3-predict.R merge|release
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  release = "C:/Users/adf44/source/r/rellib-r6",
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
set.seed(1)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(n / 10, 10), y = 0)
d$y <- frm_simulate(bf(y ~ x) + gaussian(), d,
                    newparams = list(b_Intercept = 1, b_x = 2, sigma = 1),
                    nsim = 1, seed = 1001L)[[1]]
fit <- suppressWarnings(frm(bf(y ~ x + (1 + x | g)), data = d))
cat("lost:", paste(names(frmtmb:::sdr_of(fit)$se_lost), collapse = ","), "\n")
nd <- data.frame(x = c(-3, 0, 3), g = factor(1, levels = levels(d$g)))
for (rf in list(NA, NULL)) {
  w <- character()
  p <- withCallingHandlers(
    predict(fit, newdata = nd, re_formula = rf, ndraws = 4000),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    })
  cat(sprintf("re_formula %-4s: Est.Error %s | width %s | warnings %d%s\n",
              format(rf), paste(signif(p[, "Est.Error"], 4), collapse = " "),
              paste(signif(p[, 4] - p[, 3], 4), collapse = " "), length(w),
              if (length(w)) paste0(" [", substr(w[1], 1, 90), "]") else ""))
}
f <- fitted(fit, newdata = nd, re_formula = NA)
cat("fitted() Est.Error (mean, re_formula = NA):",
    signif(f[, "Est.Error"], 4), "\n")
