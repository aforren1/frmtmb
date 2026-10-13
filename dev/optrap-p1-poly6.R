# Lane optrap, punch round 1 (m10): the review's degree-6 raw polynomial
# without a group term (probe-optrap.R poly, scale 1 and 1e5): where the
# fit stopped and how its SEs compare with lm()'s at the ML sigma.
#   Rscript dev/optrap-p1-poly6.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
for (scale in c(1, 1e5)) {
  deg <- 6
  set.seed(deg)
  d <- data.frame(x = runif(200, 1, 2) * scale, g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x / scale) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + ")))
  f <- suppressWarnings(frm(fo, data = d))
  l <- lm(fo, data = d)
  n <- nrow(d)
  p <- length(coef(l))
  se_ml <- sqrt(diag(vcov(l))) * sqrt((n - p) / n)
  se_f <- summary(f)$coefficients$mu[, "Std. Error"]
  b_f <- fixef_by_dpar(f)$mu
  cat(sprintf("scale %g units %s code %d | logLik - lm %.3g | sigma/ML %.6g\n",
              scale, deparse(f$par_units), f$opt$convergence,
              as.numeric(logLik(f)) - as.numeric(logLik(l)),
              exp(fixef_by_dpar(f)$sigma) / sqrt(mean(residuals(l)^2))))
  cat("  max |beta - lm| / SE", signif(max(abs(b_f - coef(l)) / se_ml), 3),
      "| SE / ML-SE - 1:", signif(se_f / se_ml - 1, 2), "\n")
}
