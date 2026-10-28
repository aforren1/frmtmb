# sigma under the covariance form of ar(), ma() and arma() is the
# INNOVATION sd, as in brms 2.23.0: its cholesky_cor_ar1(),
# cholesky_cor_ma1() and cholesky_cor_arma1() are the autocovariance of
# the process at unit innovation variance, scaled by sigma (row by row
# under sigma ~ x, normal_time_het_lpdf()). Through frmtmb 0.76.0 sigma
# was the marginal sd there: the same model, but every reporter, the
# standard errors, the priors and the coefficients of sigma ~ x were
# about another parameter than brms's (on the AR(1) 0.9 data of
# dev/acprior-agree.R, sigma(fit) 1.078 against brms's posterior mean
# 0.432).
#
# The reference here is written from brms's covariance, with no frmtmb
# code: at frmtmb's estimates, read as brms reads them, the dense
# multivariate normal is frmtmb's logLik and its gradient is zero.

as_data <- function(seed = 20261027, ng = 15, k = 6) {
  set.seed(seed)
  d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                  x = stats::rnorm(ng * k), x2 = stats::rnorm(ng * k))
  e <- unlist(lapply(seq_len(ng), function(i) {
    as.numeric(stats::arima.sim(list(ar = 0.7, ma = 0.3), k, sd = 0.6))
  }))
  d$y <- 1 + 0.5 * d$x + exp(0.2 * d$x2) * e
  d
}

# brms's covariance over sigma^2 for one group of k rows: the
# autocovariance of the ARMA process at unit innovation variance, from
# the psi weights of stats::ARMAtoMA()
as_brms_cov <- function(ar, ma, k) {
  psi <- c(1, stats::ARMAtoMA(ar = ar, ma = ma, lag.max = 5000))
  acv <- vapply(0:(k - 1), function(h) {
    sum(psi[1:(length(psi) - h)] * psi[(1 + h):length(psi)])
  }, 0)
  stats::toeplitz(acv)
}

# the negative log-likelihood of brms's model at natural parameters
# p = (b0, b1, log sigma intercept, [log sigma slope], ar..., ma...)
as_brms_nll <- function(d, p, np, nq, sigx) {
  b <- p[1:2]
  ls <- if (sigx) p[3:4] else p[3]
  k0 <- length(ls) + 2
  ar <- p[k0 + seq_len(np)]
  ma <- p[k0 + np + seq_len(nq)]
  mu <- b[1] + b[2] * d$x
  si <- rep_len(exp(ls[1] + if (sigx) ls[2] * d$x2 else 0), nrow(d))
  out <- 0
  for (rows in split(seq_len(nrow(d)), d$g)) {
    S <- as_brms_cov(ar, ma, length(rows)) * outer(si[rows], si[rows])
    out <- out - mvtnorm::dmvnorm(d$y[rows], mu[rows], S, log = TRUE)
  }
  out
}

test_that("sigma is brms's innovation sd: the dense reference at the
           estimates is the logLik, with zero gradient", {
  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("numDeriv")
  d <- as_data()
  cases <- list(
    list(f = y ~ x + ar(t, g, cov = TRUE), p = 1, q = 0, sigx = FALSE),
    list(f = y ~ x + ma(t, g, cov = TRUE), p = 0, q = 1, sigx = FALSE),
    list(f = y ~ x + arma(t, g, cov = TRUE), p = 1, q = 1, sigx = FALSE),
    list(f = y ~ x + ar(t, g, p = 2, cov = TRUE), p = 2, q = 0,
         sigx = FALSE),
    list(f = bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x2), p = 1, q = 0,
         sigx = TRUE))
  for (cs in cases) {
    fit <- frm(cs$f, data = d)
    ac <- fit$frame$autocor[[1L]]
    nat <- unname(autocor_natural(fit$estimates$thetaac, ac))
    p <- c(unname(fit$estimates$beta), unname(fit$estimates$betad), nat)
    nll <- function(q) as_brms_nll(d, q, cs$p, cs$q, cs$sigx)
    lab <- deparse(if (inherits(cs$f, "brmsformula")) cs$f$formula else cs$f)
    ll <- as.numeric(logLik(fit))
    expect_lt(abs(-nll(p) - ll), 1e-10 * abs(ll), label = lab)
    # the gradient of brms's model at frmtmb's estimates, against its
    # gradient at the same point with sigma read as the marginal sd, as
    # frmtmb 0.76.0 held it: the first is the optimizer's leftover
    gr <- numDeriv::grad(nll, p)
    pm <- p
    pm[3] <- p[3] + 0.5 * log(as_brms_cov(nat[seq_len(cs$p)],
                                          nat[cs$p + seq_len(cs$q)],
                                          1)[1, 1])
    gm <- numDeriv::grad(nll, pm)
    expect_lt(max(abs(gr)), 1e-3 * max(abs(gm)), label = lab)
    # every reporter reads the one parameter
    if (!cs$sigx) {
      s <- exp(fit$estimates$betad[[1L]])
      expect_identical(unname(sigma(fit)), s)
      expect_equal(summary(fit)$spec_pars["sigma", "Estimate"], s,
                   tolerance = 1e-14)
      expect_equal(unname(VarCorr(fit)$residual__$sd[1, "Estimate"]), s,
                   tolerance = 1e-14)
      expect_identical(unname(fitted(fit, dpar = "sigma")[1, "Estimate"]),
                       s)
      expect_identical(frm_linpred(fit, dpar = "sigma",
                                   type = "response")[[1L]], s)
    }
  }
})

test_that("residuals, simulate() and the row density read the marginal
           sd", {
  d <- as_data(ng = 30)
  fit <- frm(y ~ x + arma(t, g, cov = TRUE), data = d)
  ac <- fit$frame$autocor[[1L]]
  nat <- unname(autocor_natural(fit$estimates$thetaac, ac))
  # the marginal sd over sigma, from brms's covariance
  mf <- sqrt(as_brms_cov(nat[1], nat[2], 1)[1, 1])
  expect_gt(mf, 1.2)
  r <- residuals(fit)[, "Estimate"]
  expect_equal(residuals(fit, type = "pearson")[, "Estimate"],
               r / (sigma(fit) * mf), tolerance = 1e-12)
  s <- as.matrix(simulate(fit, nsim = 300, seed = 2))
  sd_sim <- stats::sd(as.vector(s - fitted(fit)[, "Estimate"]))
  # 54000 draws; the ratio to the marginal sd is 1 within a few percent,
  # where the ratio to sigma is mf
  expect_lt(abs(sd_sim / (sigma(fit) * mf) - 1), 0.1 * (mf - 1))
})

test_that("conditional_effects(method = 'posterior_predict') draws at
           brms's sigma", {
  # brms 2.23.0 drops every autocorrelation term in its
  # conditional_effects() (incl_autocor = FALSE), so its predictive band
  # is iid draws at its sigma, the innovation sd. So is frmtmb's: the
  # grid never reads the cov = TRUE term and draws at the dpar sigma.
  # Through frmtmb 0.76.0 that sigma was the marginal sd, and the band
  # was wider by 1 / sqrt(1 - ar^2), 1.82 on these data.
  set.seed(20261028)
  ng <- 20
  k <- 8
  d <- data.frame(t = rep(seq_len(k), ng), g = factor(rep(seq_len(ng),
                                                        each = k)))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + unlist(lapply(seq_len(ng), function(i) {
    as.numeric(stats::arima.sim(list(ar = 0.8), k, sd = 0.5))
  }))
  fit <- frm(y ~ x + ar(t, g, cov = TRUE), data = d)
  ac <- fit$frame$autocor[[1L]]
  arv <- unname(autocor_natural(fit$estimates$thetaac, ac))
  # the marginal sd of a row, read off the residuals (the same on both
  # versions), carried to the innovation sd by the fitted ar
  marg <- unname(residuals(fit)[1L, "Estimate"] /
                   residuals(fit, type = "pearson")[1L, "Estimate"])
  innov <- marg * sqrt(1 - arv^2)
  set.seed(1)
  ce <- conditional_effects(fit, effects = "x",
                            method = "posterior_predict",
                            ndraws = 4000)[[1L]]
  half <- mean((ce$upper__ - ce$lower__) / (2 * stats::qnorm(0.975)))
  # 4000 draws per point: the band is the innovation sd within a few
  # percent, against a factor 1 / sqrt(1 - ar^2) for the marginal sd
  expect_lt(abs(half / innov - 1), 0.1 * (1 / sqrt(1 - arv^2) - 1))
})

test_that("a fit saved with the marginal sigma is refused where it
           reads sigma", {
  d <- as_data()
  fit <- frm(y ~ x + ar(t, g, cov = TRUE), data = d)
  expect_identical(fit$frame$autocor[[1L]]$sigma_is, "innovation")
  old <- fit
  # what frmtmb 0.76.0 stored: no record of the scale, and sigma the
  # marginal sd
  old$frame$autocor[[1L]]$sigma_is <- NULL
  msg <- "frmtmb 0.76.0 or earlier"
  # every reader of sigma's value or of the coefficient vector that
  # holds it
  for (m in c("summary", "print", "sigma", "VarCorr", "confint",
              "residuals", "simulate", "predict", "cooks.distance")) {
    args <- list(old)
    if (m == "residuals") args$type <- "pearson"
    expect_error(do.call(m, args), msg, class = "frmtmb_error",
                 label = m)
  }
  expect_error(vcov(old, full = TRUE), msg, class = "frmtmb_error")
  expect_error(fitted(old, dpar = "sigma"), msg, class = "frmtmb_error")
  expect_error(frm_linpred(old, dpar = "sigma"), msg,
               class = "frmtmb_error")
  expect_error(conditional_effects(old, effects = "x", dpar = "sigma"),
               msg, class = "frmtmb_error")
  expect_error(conditional_effects(old, effects = "x",
                                   method = "posterior_predict"),
               msg, class = "frmtmb_error")
  expect_error(hypothesis(old, "sigma > 0", class = NULL), msg,
               class = "frmtmb_error")
  if (requireNamespace("insight", quietly = TRUE)) {
    expect_error(insight::get_sigma(old), msg, class = "frmtmb_error")
    expect_error(insight::get_parameters(old), msg,
                 class = "frmtmb_error")
    # and on a current fit, insight reads sigma() itself
    expect_identical(as.numeric(insight::get_sigma(fit)),
                     unname(sigma(fit)))
  }
  # a reader of the mean alone answers: the model is the same, so its
  # number is right
  expect_identical(logLik(old), logLik(fit))
  expect_identical(fitted(old), fitted(fit))
  expect_identical(residuals(old), residuals(fit))
  expect_identical(fixef(old), fixef(fit))
  expect_identical(vcov(old), vcov(fit))
  expect_identical(frm_linpred(old), frm_linpred(fit))
  # under sigma ~ x the coefficients of sigma are in fixef() and vcov()
  fx <- frm(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x2), data = d)
  oldx <- fx
  oldx$frame$autocor[[1L]]$sigma_is <- NULL
  expect_error(fixef(oldx), msg, class = "frmtmb_error")
  expect_error(vcov(oldx), msg, class = "frmtmb_error")
  expect_error(coef(oldx), msg, class = "frmtmb_error")
  # update() refits, with this version's sigma
  expect_identical(as.numeric(logLik(update(old))),
                   as.numeric(logLik(fit)))
  # the absent cases: cosy() and ar() with cov = FALSE, whose sigma was
  # and is brms's, carry no record and are not refused (the second has
  # struct "ar" too)
  for (fo in list(y ~ x + cosy(t, g), y ~ x + ar(t, g))) {
    f0 <- frm(fo, data = d)
    expect_null(f0$frame$autocor[[1L]]$sigma_is)
    expect_no_error(summary(f0))
    expect_no_error(sigma(f0))
    expect_no_error(fitted(f0, dpar = "sigma"))
  }
})
