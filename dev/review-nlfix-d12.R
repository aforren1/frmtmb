# Fix round 2, item 3: the density's cost on my data, seed 13020, n = 3e4, min of 3.
# FRM_ARM before or after; after also times my form as written, as a custom family.
source("dev/review-nlfix-common.R")
ns <- asNamespace("frmtmb")
bb2 <- if (arm == "after") frmtmb_family("bb2", dpars = c("mu", "phi"), links = list(mu = "logit", phi = "log"), type = "discrete", required_aterms = "trials",
  lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
    size <- ns$aterm_trials(aterms, "bb2"); a <- dpars$mu * dpars$phi; b <- (1 - dpars$mu) * dpars$phi
    lgamma(size + 1) - lgamma(y + 1) - lgamma(size - y + 1) + ns$lgamma_shift_diff(a, y) + ns$lgamma_shift_diff(b, size - y) - ns$lgamma_shift_diff(dpars$phi, size)
  }, init_dpars = list(mu = function(y, aterms) 0.5, phi = function(y, aterms) 5))
for (tr in c(5, 10, 24, 25, 30, 200)) {
  set.seed(13020); n <- 3e4; d <- data.frame(x = rnorm(n), ntr = as.integer(tr)); p <- plogis(0.3 + 0.5 * d$x); d$y <- rbinom(n, tr, rbeta(n, p * 8, (1 - p) * 8))
  t1 <- Inf; for (i in 1:3) t1 <- min(t1, system.time(f1 <- suppressWarnings(frm(bf(y | trials(ntr) ~ x), family = beta_binomial(), data = d)))[["elapsed"]])
  t2 <- NA; ll2 <- NA
  if (!is.null(bb2)) { t2 <- Inf; for (i in 1:2) t2 <- min(t2, system.time(f2 <- suppressWarnings(frm(bf(y | trials(ntr) ~ x), family = bb2, data = d)))[["elapsed"]]); ll2 <- as.numeric(logLik(f2)) }
  cat(sprintf("%s trials %3d: family %.2f s ll %.6f | custom form as written %.2f s ll %.6f\n", arm, tr, t1, as.numeric(logLik(f1)), t2, ll2))
}
