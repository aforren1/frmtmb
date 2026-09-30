# frmtmb.sample on the three families: frm_sample(), log_lik() against the
# reference density at the draws, posterior_epred(), posterior_predict()
# and conditional_effects(). Seed 20260929.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb.sample from", find.package("frmtmb.sample"), "\n")
say <- function(...) cat(sprintf(...), "\n", sep = "")
ref_xbeta <- function(y, mu, phi, kappa) {
  a <- mu * phi; b <- (1 - mu) * phi; d <- 1 + 2 * kappa
  ifelse(y <= 0, pbeta(kappa / d, a, b, log.p = TRUE),
         ifelse(y >= 1, pbeta((1 + kappa) / d, a, b, lower.tail = FALSE, log.p = TRUE),
                dbeta((y + kappa) / d, a, b, log = TRUE) - log(d)))
}
ref_zibb <- function(y, size, mu, phi, zi) {
  a <- mu * phi; b <- (1 - mu) * phi
  base <- lchoose(size, y) + lbeta(y + a, size - y + b) - lbeta(a, b)
  ifelse(y == 0, log(zi + (1 - zi) * exp(base)), log1p(-zi) + base)
}
ref_hc <- function(y, eta, hu, tau) {
  K1 <- length(tau)
  up <- ifelse(y >= K1 + 1, 1, plogis(tau[pmin(pmax(y, 1), K1)] - eta))
  lo <- ifelse(y <= 1, 0, plogis(tau[pmax(pmin(y - 1, K1), 1)] - eta))
  ifelse(y == 0, log(hu), log1p(-hu) + log(up - lo))
}
set.seed(20260929)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), tr = sample(4:12, n, TRUE))
mu <- plogis(0.3 + 0.6 * d$x); kap <- 0.15
zz <- rbeta(n, mu * 6, (1 - mu) * 6)
d$yx <- pmin(pmax((1 + 2 * kap) * zz - kap, 0), 1)
mu <- plogis(-0.4 + 0.5 * d$x)
d$yb <- ifelse(runif(n) < 0.25, 0L, rbinom(n, d$tr, rbeta(n, mu * 4, (1 - mu) * 4)))
u <- rlogis(n) + 0.8 * d$x
d$yh <- ifelse(runif(n) < plogis(-0.6 + 0.5 * d$z), 0L, 1L + (u > -1) + (u > 0.3) + (u > 1.5))
fits <- list(
  xbeta = frm(bf(yx ~ x), family = xbeta(), data = d),
  zibb = frm(bf(yb | trials(tr) ~ x), family = zero_inflated_beta_binomial(), data = d),
  hurdle_cumulative = frm(bf(yh ~ x, hu ~ z), family = hurdle_cumulative(), data = d))
for (nm in names(fits)) {
  dr <- withCallingHandlers(frm_sample(fits[[nm]], chains = 2, iter = 600, seed = 1, refresh = 0),
                            message = function(m) { cat("  [msg]", conditionMessage(m)); invokeRestart("muffleMessage") },
                            warning = function(w) { cat("  [warn]", conditionMessage(w), "
"); invokeRestart("muffleWarning") })
  ll <- log_lik(dr)
  S <- nrow(ll)
  pick <- unique(round(seq(1, S, length.out = 20)))
  worst <- 0
  for (s in pick) {
    f <- frmtmb.sample:::draws_fit_at(dr, s, frmtmb.sample:::draws_par_index(dr$fit))
    dp <- frmtmb:::eval_dpars(f)[[1]]
    ref <- switch(nm,
      xbeta = ref_xbeta(d$yx, dp$mu, dp$phi, dp$kappa),
      zibb = ref_zibb(d$yb, d$tr, dp$mu, dp$phi, dp$zi),
      hurdle_cumulative = ref_hc(d$yh, dp$mu, dp$hu,
                                 frmtmb:::ord_tau_from_raw(f$estimates$tau_raw, TRUE)))
    worst <- max(worst, max(abs(ll[s, ] - ref) / pmax(1, abs(ref))))
  }
  say("%-18s %d draws; log_lik vs reference at %d draws: max rel diff %.2e",
      nm, S, length(pick), worst)
  ep <- posterior_epred(dr, ndraws = 50)
  say("  posterior_epred dim %s", paste(dim(ep), collapse = " x "))
  pp <- posterior_predict(dr, ndraws = 50)
  say("  posterior_predict dim %s, range %s", paste(dim(pp), collapse = " x "),
      paste(range(pp), collapse = " .. "))
  ce <- tryCatch({conditional_effects(dr, effects = "x", ndraws = 50); "ok"},
                 error = function(e) conditionMessage(e))
  say("  conditional_effects: %s", ce)
  rh <- tryCatch(max(posterior::summarise_draws(posterior::as_draws_array(dr), "rhat")$rhat, na.rm = TRUE), error = function(e) NA)
  say("  max rhat %.3f", rh)
  lo <- tryCatch({l <- loo(dr); sprintf("elpd_loo %.1f", l$estimates["elpd_loo", 1])},
                 error = function(e) conditionMessage(e))
  say("  loo: %s", lo)
}
