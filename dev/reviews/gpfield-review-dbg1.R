# Review of lane gpfield: debug of gpfield-review-kernels.R's NA rows, the
# matern32 derivative, and the warnings seen.
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
suppressPackageStartupMessages(library(brms))
ns <- asNamespace("frmtmb"); bns <- asNamespace("brms")
mk <- function(seed, n) {
  set.seed(seed)
  x1 <- round(runif(n, 0, 10), 2); x2 <- round(runif(n, 0, 5), 2)
  g <- factor(sample(c("a", "b", "c"), n, TRUE))
  y <- sin(x1) * (1 + (g == "b")) + 0.3 * x2 + 0.2 * rnorm(n)
  data.frame(x1, x2, g, y)
}
d <- mk(9404, 45L)
fml <- y ~ gp(x1, x2, cov = "exponential")
sdat <- brms::standata(brms::bf(fml), data = d, family = gaussian())
print(names(sdat))
q <- gf_fit(frmtmb::bf(fml), data = d, family = gaussian())
f <- q$fit
print(f$estimates$beta); print(f$estimates$betad); print(f$estimates$theta)
d <- mk(9409, 45L)
q <- gf_fit(frmtmb::bf(y ~ gp(x1, x2, k = 7, c = c(1.5, 2), iso = FALSE,
                           cov = "matern52")), data = d, family = gaussian())
print(q$warn)
# matern32 derivative, by case 1, with a range of steps
d <- mk(9401, 45L)
q <- gf_fit(frmtmb::bf(y ~ g + gp(x1, by = g, cov = "matern32")), data = d,
            family = gaussian())
f <- q$fit
for (bk in Filter(function(b) b$covstruct == "gp", f$frame$re_blocks)) {
  th <- f$estimates$theta[bk$theta_idx]
  z <- f$estimates$b[bk$b_idx]
  L0 <- ns$gp_factor_ordered(th, bk)
  an <- exp(th[1]) * (L0 %*% ns$gp_field_dtheta(z, th, bk))[, 2]
  K <- function(t) { k <- ns$gp_corr(t, bk); k }
  # dK check directly: analytic dK from the kernel against central diff
  r <- sqrt(Reduce(`+`, bk$aux_D2))
  dKa <- ns$gp_kern_dlog(r / exp(th[2]), "matern32")
  for (h in c(1e-3, 1e-4, 1e-5, 1e-6)) {
    tp <- th; tp[2] <- th[2] + h; tm <- th; tm[2] <- th[2] - h
    dKn <- (K(tp) - K(tm)) / (2 * h)
    fd <- exp(th[1]) * (ns$gp_factor_ordered(tp, bk) %*% z -
                          ns$gp_factor_ordered(tm, bk) %*% z) / (2 * h)
    cat(sprintf("%s h %.0e | dK rel err %.2e | dfield rel err %.2e | kappa %.2e ls %.3f npos %d\n",
                bk$group_name, h, max(abs(dKn - dKa)) / max(abs(dKa)),
                max(abs(fd - an)) / max(abs(fd)),
                kappa(K(th)), exp(th[2]), nrow(r)))
  }
}
