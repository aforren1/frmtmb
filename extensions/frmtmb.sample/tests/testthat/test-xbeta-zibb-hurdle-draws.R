# The draws surface on core's xbeta(), zero_inflated_beta_binomial() and
# hurdle_cumulative(). log_lik() is checked against densities written
# here from stats::pbeta(), stats::dbeta(), lbeta() and plogis(), not
# from the package code; dev/fams2-sample-smoke.R prints the same
# comparisons.

f2_ref_xbeta <- function(y, mu, phi, kappa) {
  a <- mu * phi
  b <- (1 - mu) * phi
  d <- 1 + 2 * kappa
  ifelse(y <= 0, stats::pbeta(kappa / d, a, b, log.p = TRUE),
         ifelse(y >= 1,
                stats::pbeta((1 + kappa) / d, a, b, lower.tail = FALSE,
                             log.p = TRUE),
                stats::dbeta((y + kappa) / d, a, b, log = TRUE) - log(d)))
}

f2_ref_zibb <- function(y, size, mu, phi, zi) {
  a <- mu * phi
  b <- (1 - mu) * phi
  base <- lchoose(size, y) + lbeta(y + a, size - y + b) - lbeta(a, b)
  ifelse(y == 0, log(zi + (1 - zi) * exp(base)), log1p(-zi) + base)
}

f2_ref_hc <- function(y, eta, hu, tau) {
  K1 <- length(tau)
  a <- c(tau, Inf)[pmin(pmax(y, 1), K1 + 1)] - eta
  b <- c(-Inf, tau)[pmin(pmax(y, 1), K1 + 1)] - eta
  p <- ifelse(a + b > 0,
              stats::plogis(b, lower.tail = FALSE) -
                stats::plogis(a, lower.tail = FALSE),
              stats::plogis(a) - stats::plogis(b))
  ifelse(y == 0, log(hu), log1p(-hu) + log(p))
}

f2_data <- function(seed, n = 150) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n), tr = sample(4:12, n, TRUE))
  mu <- stats::plogis(0.3 + 0.6 * d$x)
  zz <- stats::rbeta(n, mu * 6, (1 - mu) * 6)
  d$yx <- pmin(pmax(1.3 * zz - 0.15, 0), 1)
  mu <- stats::plogis(-0.4 + 0.5 * d$x)
  d$yb <- ifelse(runif(n) < 0.25, 0L,
                 stats::rbinom(n, d$tr, stats::rbeta(n, mu * 4,
                                                     (1 - mu) * 4)))
  u <- stats::rlogis(n) + 0.8 * d$x
  d$yh <- ifelse(runif(n) < stats::plogis(-0.6 + 0.5 * d$z), 0L,
                 1L + (u > -1) + (u > 0.3) + (u > 1.5))
  d
}

f2_draws <- local({
  cache <- list()
  function(which) {
    skip_sampler()
    if (is.null(cache[[which]])) {
      d <- f2_data(61)
      fit <- switch(which,
        xbeta = frm(bf(yx ~ x), family = xbeta(), data = d),
        zibb = frm(bf(yb | trials(tr) ~ x),
                   family = zero_inflated_beta_binomial(), data = d),
        hurdle = frm(bf(yh ~ x, hu ~ z), family = hurdle_cumulative(),
                     data = d))
      msg <- character()
      ds <- withCallingHandlers(
        suppressWarnings(frm_sample(fit, chains = 1, iter = 300,
                                    refresh = 0, seed = 3)),
        message = function(m) {
          msg <<- c(msg, conditionMessage(m))
          invokeRestart("muffleMessage")
        })
      cache[[which]] <<- list(d = d, fit = fit, ds = ds, msg = msg)
    }
    cache[[which]]
  }
})

test_that("log_lik() is the reference density at the draws", {
  for (which in c("xbeta", "zibb", "hurdle")) {
    cs <- f2_draws(which)
    ll <- log_lik(cs$ds)
    idx <- frmtmb.sample:::draws_par_index(cs$ds$fit)
    for (s in c(1L, nrow(ll))) {
      f <- frmtmb.sample:::draws_fit_at(cs$ds, s, idx)
      dp <- frmtmb:::eval_dpars(f)[[1]]
      ref <- switch(which,
        xbeta = f2_ref_xbeta(cs$d$yx, dp$mu, dp$phi, dp$kappa),
        zibb = f2_ref_zibb(cs$d$yb, cs$d$tr, dp$mu, dp$phi, dp$zi),
        hurdle = f2_ref_hc(cs$d$yh, dp$mu, dp$hu,
                           frmtmb:::ord_tau_from_raw(f$estimates$tau_raw,
                                                     TRUE)))
      expect_lt(max(abs(ll[s, ] - ref) / pmax(1, abs(ref))),
                1e3 * .Machine$double.eps, label = which)
    }
  }
})

test_that("a hurdle ordinal's expected category is scored 0..K", {
  cs <- f2_draws("hurdle")
  pp <- posterior_predict(cs$ds, ndraws = 20)
  expect_true(all(pp %in% 0:4))
  ep <- posterior_epred(cs$ds, ndraws = 20)
  expect_identical(dim(ep)[3], 5L)
  cd <- allow_warnings(
    conditional_effects(cs$ds, effects = "x", resolution = 4,
                        categorical = FALSE, robust = FALSE),
    "Predictions are treated as continuous")
  pc <- conditional_effects(cs$ds, effects = "x", resolution = 4,
                            robust = FALSE)[["x:cats__"]]
  ncat <- nlevels(pc$cats__)
  expect_identical(ncat, 5L)
  ngrid <- nrow(pc) / ncat
  # the codes, 0 for the hurdle, and not the column positions
  by_hand <- rowSums(vapply(seq_len(ncat), function(k) {
    (k - 1) * pc$estimate__[(k - 1L) * ngrid + seq_len(ngrid)]
  }, numeric(ngrid)))
  expect_lt(max(abs(cd$x$estimate__ - by_hand)), 1e-10)
})

test_that("the default-prior disclosure names xbeta's kappa", {
  cs <- f2_draws("xbeta")
  expect_true(any(grepl("no defaults for phi, kappa", cs$msg,
                        fixed = TRUE)))
  cz <- f2_draws("zibb")
  expect_true(any(grepl("no defaults for phi", cz$msg, fixed = TRUE)))
})
