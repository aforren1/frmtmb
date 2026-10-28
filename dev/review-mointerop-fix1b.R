## Reviewer: non-gap sandwiches on rellib-r14 (base), compared with the
## lane's saved values from dev/review-mointerop-fix1.R.
##   Rscript dev/review-mointerop-fix1b.R base
source("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-common.R")
set.seed(31)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(sample(1:30, n, TRUE)),
                w = runif(n, 0.5, 2), off = log(runif(n, 1, 3)))
d$ypos <- rgamma(n, shape = 2, rate = 2 / exp(0.5 + 0.3 * d$x))
d$ybin <- rbinom(n, 1, pnorm(0.2 + 0.6 * d$x))
d$cnt <- rpois(n, exp(0.2 + 0.4 * d$x + d$off))
d$y <- 1 + 0.5 * d$x + rnorm(n) * (1 + abs(d$x))
f0 <- frm(y ~ x, data = d)
fp <- frm(cnt ~ x + offset(off), data = d, family = poisson())
dm <- data.frame(k = sample(1:4, n, TRUE), z = rnorm(n),
                 g = factor(sample(1:30, n, TRUE)))
dm$y <- 1 + c(0, 0.3, 0.5, 0.9)[dm$k] + 0.4 * dm$z + rnorm(n)
fm <- frm(y ~ mo(k) + z, data = dm)
set.seed(33)
dn <- data.frame(x = rnorm(n), g = factor(rep(1:30, each = 10)))
dn$y <- 1 + 0.5 * dn$x + rnorm(30)[dn$g] + rnorm(n)
fn <- frm(y ~ x + (1 | g), data = dn)
base <- list(
  f0_hc0 = vcov_cluster(f0, seq_len(n), "CR0", full = TRUE),
  f0_g = vcov_cluster(f0, ~ g, "CR1", full = TRUE),
  fp_g = vcov_cluster(fp, ~ g, "CR1", full = TRUE),
  fm_g = vcov_cluster(fm, ~ g, "CR1", full = TRUE),
  fn_g = vcov_cluster(fn, ~ g, "CR1", full = TRUE),
  f0_me = vcov_cluster(f0, seq_len(n), "CR0"),
  fp_me = vcov_cluster(fp, ~ g, "CR1"),
  fn_me = vcov_cluster(fn, ~ g, "CR1"))
lane <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/review-mointerop-out/fix1-nongap-lane.rds")
for (k in names(base)) {
  a <- unclass(base[[k]]); b <- unclass(lane[[k]])
  attributes(a) <- attributes(a)["dim"]; attributes(b) <- attributes(b)["dim"]
  cat(sprintf("%-7s identical %-5s max rel diff %.3g\n", k,
              identical(a, b), max(abs(a - b)) / max(abs(a))))
}
