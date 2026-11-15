# Lane gpfield: predict()'s spread at x = 5.5 on an approximate field
# (test-gp-predict-draw.R's data), on one arm. Args: arm
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs(a[1])
set.seed(9200)
n <- 30
d <- data.frame(x = sort(stats::runif(n, 0, 10)))
d$y <- 2 * sin(d$x) + stats::rnorm(n, 0, 0.2)
q <- gf_fit(bf(y ~ gp(x, k = 15)), data = d, family = gaussian())
fit <- q$fit
cat("warn", q$warn, "\n")
cat("theta", fit$estimates$theta, "sigma", exp(fit$estimates$betad), "\n")
print(sqrt(diag(vcov(fit, full = TRUE))))
for (pe in c(TRUE, FALSE)) {
  set.seed(2)
  p <- predict(fit, newdata = data.frame(x = c(5.5, 60)), summary = FALSE,
               ndraws = 2000, propagate_error = pe)
  cat("propagate_error", pe, "var", apply(p, 2, var), "quantiles",
      quantile(p[, 1], c(0.001, 0.5, 0.999)), "\n")
}
f2 <- frm(bf(y ~ gp(x)), data = d, family = gaussian())
set.seed(1)
p <- predict(f2, newdata = data.frame(x = c(60, d$x[5])), summary = FALSE,
             ndraws = 2000)
cat("exact: sd far", sd(p[, 1]), "sd at x5", sd(p[, 2]), "sdgp",
    exp(f2$estimates$theta[1]), "sigma", exp(f2$estimates$betad[1]), "\n")
