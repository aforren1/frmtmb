# Reviewer of lane setier: copy of dev/nanse-rev-grby.R with lane = wt-setier-lib over
# rellib-r6, base = rellib-r6 (merge arm = rellib-r6 alone, i.e. base).
# Reviewer: frmtmb.sample's gr(g, by = f) fixture, where the lane names
# x, theta_2 and theta_3 "concave" and base left 2 of 9 SEs non-finite.
# Was base's finite SE of x meaningful? References: lme4 with the same
# covariance per by-level (dummy coding), and the base sdreport.
#   Rscript dev/nanse-rev-grby.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
set.seed(11)
dd <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, 10)))
dd$f <- factor(ifelse(as.integer(dd$g) <= 8, "a", "b"))
u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
dd$y <- stats::rnorm(160, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
fit <- suppressWarnings(frm(bf(y ~ x + (1 + x | gr(g, by = f))),
                            family = gaussian(), data = dd))
nm <- ns$outer_par_names(fit)
cat("logLik", format(logLik(fit), digits = 12), "code",
    fit$opt$convergence, "\n")
cat("par:", paste(sprintf("%s=%.4g", nm, fit$opt$par), collapse = " "), "\n")
V <- ns$sdr_of(fit)$cov.fixed
cat("reported SE:", paste(sprintf("%s=%.4g", nm,
                                  suppressWarnings(sqrt(diag(V)))),
                          collapse = " "), "\n")
H <- stats::optimHess(fit$opt$par, fit$obj$fn, fit$obj$gr)
H <- (H + t(H)) / 2
D <- sqrt(abs(diag(H)))
e <- eigen(H / outer(D, D), symmetric = TRUE)
cat("unit-diag eigenvalues:", paste(signif(e$values, 3), collapse = " "), "\n")
for (k in which(e$values < 1e-3 * max(e$values))) {
  cat(sprintf("  eigenvalue %.3g loads:", e$values[k]),
      paste(sprintf("%s=%.2f", nm, e$vectors[, k])[abs(e$vectors[, k]) >
                                                     0.05],
            collapse = " "), "\n")
}
if (requireNamespace("lme4", quietly = TRUE)) {
  dd$fa <- as.numeric(dd$f == "a")
  dd$fb <- as.numeric(dd$f == "b")
  m <- suppressWarnings(suppressMessages(lme4::lmer(
    y ~ x + (0 + fa + fa:x | g) + (0 + fb + fb:x | g), data = dd,
    REML = FALSE)))
  cat("lme4 logLik", format(logLik(m), digits = 12), " SE",
      paste(sprintf("%s=%.4g", names(lme4::fixef(m)),
                    sqrt(diag(as.matrix(vcov(m))))), collapse = " "),
      " singular", lme4::isSingular(m), "\n")
  print(lme4::VarCorr(m))
}
