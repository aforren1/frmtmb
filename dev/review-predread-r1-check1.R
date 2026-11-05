# Short check, item 1: separate Monte Carlo error from delta-method error on the rows
# left on the M1 to M3 table. Seed 7971, the data of review-predread-mefam.R. For each
# family, at x = 1, sdx = 1.2: the closed form or a 2e6-draw integral over the latent at
# the fit's estimate, 50 calls of brms's posterior_epred at fixed values (200000 draws),
# and frmtmb's Estimate, Est.Error with the coefficient part shown.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7971)
n <- 400
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4)
d$xo <- d$lat + rnorm(n, 0, 0.4)
eta <- 0.4 + 0.7 * d$lat
d$cnt <- rpois(n, exp(eta))
d$nb <- rnbinom(n, mu = exp(eta), size = 3)
d$pos <- exp(eta + rnorm(n, 0, 0.5))
d$gam <- rgamma(n, 4, 4 / exp(eta))
d$frac <- rbeta(n, plogis(eta) * 8, (1 - plogis(eta)) * 8)
d$bin <- rbinom(n, 1, plogis(eta))
d$zip <- ifelse(runif(n) < 0.3, 0L, d$cnt)
d$hp <- ifelse(d$cnt == 0, 1L, d$cnt)
d$hp[runif(n) < 0.3] <- 0L
d$ord <- as.integer(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE))
d$cat <- as.integer(cut(eta + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
nd <- data.frame(xo = 1, sx = 1.2)
fmt <- function(v) paste(signif(as.numeric(v), 5), collapse = " ")
rep_epred <- function(fx, K = 50L, cat = NULL) {
  out <- vector("list", K)
  for (i in seq_len(K)) {
    a <- posterior_epred(fx, newdata = nd)
    out[[i]] <- if (length(dim(a)) == 3L) a[, 1, cat] else a[, 1]
  }
  unlist(out)
}
mcse_sd <- function(v) {
  m4 <- mean((v - mean(v))^4)
  sqrt((m4 - var(v)^2) / length(v)) / (2 * sd(v))
}
mono <- function(nm, sp, fam, tr) {
  cat("\n===", nm, "\n")
  pr <- rev_pair(sp, fam, d, ndraws = 4000)
  fit <- pr$fit
  b <- fit$estimates$beta
  print(signif(unlist(fit$estimates[c("beta", "betad")]), 5))
  bi <- b[1]
  bs <- b[grep("mexosx", names(b))]
  m <- bi + bs * nd$xo
  s <- abs(bs) * nd$sx
  set.seed(11)
  cf <- tr(fit, m, s)
  v <- rep_epred(pr$fx)
  fv <- fitted(fit, newdata = nd)
  lb <- frm_lp_basis(fit, newdata = nd)
  avar <- rowSums((lb$A %*% lb$V) * lb$A)
  cat("closed form mean, sd:", fmt(cf), "\n")
  cat("brms 200000 draws mean, sd:", fmt(c(mean(v), sd(v))), " MCSE of mean, sd:", fmt(c(sd(v) / sqrt(length(v)), mcse_sd(v))), "\n")
  first <- v[seq_len(4000)]
  cat("brms first 4000 draws sd:", signif(sd(first), 4), " MCSE:", signif(mcse_sd(first), 3), "\n")
  sds4k <- vapply(split(v, rep(seq_len(50), each = 4000)), sd, 0)
  cat("sd over 50 sets of 4000: min", signif(min(sds4k), 4), " median", signif(median(sds4k), 4), " max", signif(max(sds4k), 4), "\n")
  cat("frmtmb Estimate, Est.Error:", fmt(fv[, c("Estimate", "Est.Error")]), " Q2.5/Q97.5:", fmt(fv[, c("Q2.5", "Q97.5")]), "\n")
  cat("coef var on eta (AVA'):", signif(avar, 4), " eta draw var:", signif(s^2, 4), "\n")
  cat("brms quantiles:", fmt(quantile(v, c(0.025, 0.975))), "\n")
  invisible(list(fit = fit, fx = pr$fx))
}
ln <- function(m, s) c(exp(m + s^2 / 2), exp(m + s^2 / 2) * sqrt(exp(s^2) - 1))
mono("gamma", "bf(gam ~ me(xo, sx))", "Gamma('log')", function(fit, m, s) ln(m, s))
mono("poisson", "bf(cnt ~ me(xo, sx))", "poisson()", function(fit, m, s) ln(m, s))
mono("zero_inflated_poisson", "bf(zip ~ me(xo, sx))", "zero_inflated_poisson()", function(fit, m, s) {
  zi <- plogis(fit$estimates$betad[grep("zi", names(fit$estimates$betad))[1L]])
  cat("zi:", signif(zi, 5), "\n")
  (1 - zi) * ln(m, s)
})
cat("\n=== categorical: does the mean of category 1 move with sharing the latent?\n")
pr <- rev_pair("bf(cat ~ me(xo, sx))", "categorical()", d, ndraws = 4000)
fit <- pr$fit
b <- fit$estimates$beta
print(signif(b, 5))
a2 <- b[grep("mu2.*Intercept", names(b))]
a3 <- b[grep("mu3.*Intercept", names(b))]
b2 <- b[grep("mu2.*me", names(b))]
b3 <- b[grep("mu3.*me", names(b))]
M <- 2e6
set.seed(12)
z1 <- rnorm(M)
z2 <- rnorm(M)
pcat <- function(zz2, zz3) {
  e2 <- a2 + b2 * (nd$xo + nd$sx * zz2)
  e3 <- a3 + b3 * (nd$xo + nd$sx * zz3)
  den <- 1 + exp(e2) + exp(e3)
  cbind(1 / den, exp(e2) / den, exp(e3) / den)
}
sh <- pcat(z1, z1)
ind <- pcat(z1, z2)
cat("one latent read by mu2, mu3:   mean", fmt(colMeans(sh)), " sd", fmt(apply(sh, 2, sd)), "\n")
cat("one latent per predictor:      mean", fmt(colMeans(ind)), " sd", fmt(apply(ind, 2, sd)), "\n")
cat("plug-in at x:                  ", fmt(pcat(0, 0)), "\n")
for (k in 1:3) {
  v <- rep_epred(pr$fx, K = 25L, cat = k)
  cat("brms category", k, "100000 draws: mean", signif(mean(v), 5), " sd", signif(sd(v), 5), "\n")
}
fv <- fitted(fit, newdata = nd)
cat("frmtmb fitted mean:", fmt(fv[, "Estimate", ]), " Est.Error:", fmt(fv[, "Est.Error", ]), "\n")
cat("frmtmb fitted Q2.5:", fmt(fv[, "Q2.5", ]), " Q97.5:", fmt(fv[, "Q97.5", ]), "\n")
ds <- rev_draws_at(fit, 20000L)
set.seed(13)
e <- frmtmb.sample::posterior_epred(ds, newdata = nd)
cat("frmtmb draws posterior_epred mean:", fmt(colMeans(e[, 1, ])), " sd:", fmt(apply(e[, 1, ], 2, sd)), "\n")
cat("frmtmb draws Q2.5:", fmt(apply(e[, 1, ], 2, quantile, 0.025)), " Q97.5:", fmt(apply(e[, 1, ], 2, quantile, 0.975)), "\n")
cat("\n=== cumulative: the clamped Wald interval against the quantile, x = 1, sdx = 1.2\n")
pr2 <- rev_pair("bf(ord ~ me(xo, sx))", "cumulative()", d, ndraws = 4000)
fv <- fitted(pr2$fit, newdata = nd)
print(signif(fv[1, , ], 4))
for (k in 1:4) {
  v <- rep_epred(pr2$fx, K = 25L, cat = k)
  cat("brms category", k, "100000 draws: mean", signif(mean(v), 4), " sd", signif(sd(v), 4), " q", fmt(quantile(v, c(0.025, 0.975))), "\n")
}
ceo <- conditional_effects(pr2$fit, "xo", resolution = 2, int_conditions = list(xo = 1), conditions = data.frame(sx = 1.2))[[1]]
print(as.data.frame(ceo)[, intersect(c("cats__", "estimate__", "lower__", "upper__"), names(ceo))], digits = 4)
p <- fv[1, "Estimate", ]
se <- fv[1, "Est.Error", ]
lp <- qlogis(p)
sel <- se / (p * (1 - p))
cat("logit-scale band of the same Estimate, Est.Error, category by category:\n")
print(signif(rbind(lo = plogis(lp - 1.96 * sel), hi = plogis(lp + 1.96 * sel)), 4))
