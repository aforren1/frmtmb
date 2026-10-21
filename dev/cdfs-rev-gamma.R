# Reviewer: where gamma_cdf_log() is NaN at large shapes, and whether a
# censored Gamma fit can reach it.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
for (a in c(1e4, 3e4, 5e4, 1e5, 1e6)) {
  sd <- sqrt(a)
  x <- a + sd * c(-8, -6, -5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5, 6, 8)
  p <- ns$gamma_cdf_log(x / a, 1, a)
  rF <- stats::pgamma(x, a, log.p = TRUE)
  rS <- stats::pgamma(x, a, lower.tail = FALSE, log.p = TRUE)
  cat(sprintf("shape %g\n", a))
  print(data.frame(k_sd = (x - a) / sd, lF = p$lF, ref_lF = rF, lS = p$lS,
                   ref_lS = rS), digits = 10)
}
cat("== the pieces at a = 1e5, x one sd below the mode\n")
a <- 1e5
x <- a - sqrt(a)
g <- ns$igamma_gl_tails(x, a)
print(g)
cat("lgamma_binet(1e5):", ns$lgamma_binet(1e5), "\n")
cat("lgamma_binet(1e4):", ns$lgamma_binet(1e4), "\n")
cat("lgamma_binet(3e4):", ns$lgamma_binet(3e4), "\n")
print(body(ns$lgamma_binet))
cat("== a censored Gamma fit with shape near 1e5 (cv 0.3%)\n")
set.seed(3)
d <- data.frame(x = rnorm(200))
d$y <- rgamma(200, shape = 1e5, rate = 1e5 / exp(1 + 0.1 * d$x))
d$cc <- as.integer(d$y > quantile(d$y, 0.8))
d$y[d$cc == 1] <- quantile(d$y, 0.8)
r <- tryCatch({
  fit <- frm(y | cens(cc) ~ x, family = Gamma(link = "log"), data = d)
  c(conv = fit$opt$convergence, shape = exp(coef(fit)[["shape_Intercept"]] %||% NA),
    ll = as.numeric(logLik(fit)))
}, error = function(e) conditionMessage(e))
print(r)
fitu <- tryCatch(frm(y ~ x, family = Gamma(link = "log"), data = d),
                 error = function(e) conditionMessage(e))
if (!is.character(fitu)) print(summary(fitu)) else print(fitu)
