# Item 5: prior = "flat" on a smooth. brms's and frm_sample()'s 90th percentile
# of sds differ (19.7 and 11.7, review-refuse-misc5.R). Which is right? An exact
# reference: for this gaussian model the coefficients integrate out, so the
# posterior of (sds, sigma) is on a 2-D grid. Flat on sds and sigma, as brms's
# "" rows. Then 4 chains x 6000 draws on each sampler, with n_eff. Seed 10701.
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
set.seed(10701)
n <- 150
d <- data.frame(x = runif(n), z = runif(n)); d$y <- sin(6 * d$x) + rnorm(n, 0, 0.4)
fb <- brms::bf(y ~ s(x, k = 5))
pri <- E$brms_flat_prior(fb, data = d, family = stats::gaussian())
code <- brms::make_stancode(fb, data = d, prior = pri)
sdat <- E$brms_standata(fb, data = d, prior = pri)
X <- cbind(1, sdat$Xs); Z <- sdat$Zs_1_1; y <- sdat$Y
pr("Z dim %s; X dim %s", paste(dim(Z), collapse = "x"), paste(dim(X), collapse = "x"))
ZZ <- tcrossprod(Z)
lgrid_s <- seq(log(0.2), log(5000), length.out = 241)
lgrid_e <- seq(log(0.25), log(0.8), length.out = 81)
lp <- matrix(NA_real_, length(lgrid_s), length(lgrid_e))
eZ <- eigen(ZZ, symmetric = TRUE)
U <- eZ$vectors; ev <- pmax(eZ$values, 0)
Uy <- crossprod(U, y); UX <- crossprod(U, X)
for (i in seq_along(lgrid_s)) for (j in seq_along(lgrid_e)) {
  v <- exp(2 * lgrid_e[j]) + exp(2 * lgrid_s[i]) * ev
  A <- crossprod(UX / sqrt(v))
  b <- solve(A, crossprod(UX / v, Uy))
  r <- Uy - UX %*% b
  lp[i, j] <- -0.5 * sum(log(v)) - 0.5 * as.numeric(determinant(A)$modulus) -
    0.5 * sum(r^2 / v) + lgrid_s[i] + lgrid_e[j]
}
w <- exp(lp - max(lp))
ms <- rowSums(w); ms <- ms / sum(ms)
cdf <- cumsum(ms)
qs <- c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
ex <- vapply(qs, function(p) exp(approx(cdf, lgrid_s, p)$y), 0)
pr("exact   sds quantiles %s; mass above the grid's top %.1e; sigma mean %.4f", paste(round(ex, 3), collapse = " "), ms[length(ms)],
   sum(colSums(w) / sum(w) * exp(lgrid_e)))
sb <- suppressWarnings(rstan::sampling(E$brms_stan_model(code), data = sdat, chains = 4, iter = 4000, refresh = 0, seed = 11, cores = 1))
eb <- rstan::extract(sb)
pr("brms    sds quantiles %s; n_eff %.0f, Rhat %.3f; sigma %.4f", paste(round(quantile(eb$sds_1, qs), 3), collapse = " "),
   rstan::summary(sb)$summary["sds_1[1]", "n_eff"], rstan::summary(sb)$summary["sds_1[1]", "Rhat"], mean(eb$sigma))
for (rp in c(TRUE, FALSE)) {
  sf <- suppressMessages(suppressWarnings(frm_sample(bf(y ~ s(x, k = 5)), data = d, prior = "flat", chains = 4, iter = 4000, seed = 11, refresh = 0, cores = 1, reparameterize = rp)))
  dm <- posterior::as_draws_df(sf)
  v <- dm[["sds_sx_1"]]
  pr("frm_sample reparameterize=%s sds quantiles %s; ess_tail %.0f rhat %.3f; sigma %.4f", rp, paste(round(quantile(v, qs), 3), collapse = " "),
     posterior::ess_tail(posterior::extract_variable_matrix(sf, "sds_sx_1")), posterior::rhat(posterior::extract_variable_matrix(sf, "sds_sx_1")), mean(dm[["sigma"]]))
}
# what the table and the code say the flat sds is
print(as.data.frame(frmtmb::default_prior(bf(y ~ s(x, k = 5)), data = d, route = "sample"))[, c("prior", "class", "coef", "lb")])
