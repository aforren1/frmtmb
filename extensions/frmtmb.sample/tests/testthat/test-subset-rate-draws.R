# The posterior methods on models with brms's subset() and rate():
# one response at a time under subset(), as brms asks, and the exposure
# of rate() read from newdata. The densities are checked against the
# closed form at a draw, as test-loo.R checks them.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

sr_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(27)
      n <- 60
      d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                      s1 = rep(c(TRUE, FALSE), n / 2),
                      s2 = c(rep(TRUE, 40), rep(FALSE, 20)),
                      time = stats::runif(n, 0.5, 3))
      d$y1 <- 1 + d$x + stats::rnorm(n)
      d$y2 <- stats::rpois(n, exp(0.3 - 0.4 * d$z) * d$time)
      d$y2[!d$s2] <- NA
      d$z[!d$s2] <- NA
      fit <- frm(bf(y1 | subset(s1) ~ x) + gaussian() +
                   bf(y2 | subset(s2) + rate(time) ~ z) + poisson(),
                 data = d)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 3)))
      cache <<- list(d = d, fit = fit, ds = ds)
    }
    cache
  }
})

test_that("log_lik() of a subset model is one response's own rows", {
  cs <- sr_case()
  msg <- "argument 'resp' must be a single variable name"
  expect_error(log_lik(cs$ds), msg)
  expect_error(loo(cs$ds), msg)
  # brms's nobs(): the data's rows, and one response's with resp
  expect_identical(nobs(cs$ds), nrow(cs$d))
  expect_identical(nobs(cs$ds, resp = "y1"), sum(cs$d$s1))
  ll1 <- log_lik(cs$ds, resp = "y1")
  ll2 <- log_lik(cs$ds, resp = "y2")
  expect_identical(dim(ll1), c(nrow(cs$ds$draws), sum(cs$d$s1)))
  expect_identical(dim(ll2), c(nrow(cs$ds$draws), sum(cs$d$s2)))
  idx <- frmtmb.sample:::draws_par_index(cs$ds$fit)
  d1 <- cs$d[cs$d$s1, ]
  d2 <- cs$d[cs$d$s2, ]
  bn <- names(cs$fit$frame$par_template$beta)
  for (i in c(1L, 50L)) {
    sh <- frmtmb.sample:::draws_fit_at(cs$ds, i, idx)
    b <- sh$estimates$beta
    mu1 <- b[bn == "y1_(Intercept)"] + b[bn == "y1_x"] * d1$x
    sig <- exp(sh$estimates$betad[1])
    ref1 <- stats::dnorm(d1$y1, mu1, sig, log = TRUE)
    expect_lt(max(abs(ll1[i, ] - ref1)), 1e3 * .Machine$double.eps *
                max(abs(ref1)))
    # rate(time): the poisson mean is exp(eta) * time
    mu2 <- exp(b[bn == "y2_(Intercept)"] + b[bn == "y2_z"] * d2$z) * d2$time
    ref2 <- stats::dpois(d2$y2, mu2, log = TRUE)
    expect_lt(max(abs(ll2[i, ] - ref2)), 1e3 * .Machine$double.eps *
                max(abs(ref2)))
  }
})

test_that("the posterior predictions take one response and its rows", {
  cs <- sr_case()
  msg <- "argument 'resp' must be a single variable name"
  expect_error(posterior_epred(cs$ds, ndraws = 5), msg)
  expect_error(posterior_predict(cs$ds, ndraws = 5), msg)
  expect_error(posterior_linpred(cs$ds, ndraws = 5), msg)
  expect_identical(ncol(posterior_epred(cs$ds, resp = "y1", ndraws = 5)),
                   sum(cs$d$s1))
  nd <- cs$d[1:10, ]
  expect_identical(ncol(posterior_predict(cs$ds, resp = "y1",
                                          newdata = nd, ndraws = 5)),
                   sum(nd$s1))
  # the exposure is newdata's own: the expected count over the mu dpar
  # is newdata's time, row by row
  nd$time <- seq(1, 10)
  ep <- posterior_epred(cs$ds, resp = "y2", newdata = nd, ndraws = 5)
  lp <- posterior_linpred(cs$ds, resp = "y2", newdata = nd, ndraws = 5,
                          transform = TRUE)
  keep <- nd$time[nd$s2]
  r <- sweep(ep / lp, 2, keep)
  expect_lt(max(abs(r)), 64 * .Machine$double.eps * max(keep))
  pp <- posterior_predict(cs$ds, resp = "y2", newdata = nd, ndraws = 50)
  expect_identical(ncol(pp), sum(nd$s2))
  # the draws are made around mu * time, not mu: at an exposure of 1e4
  # a draw without it would be four orders of magnitude off
  nd2 <- nd[nd$s2, ][1:2, ]
  nd2$time <- c(1, 1e4)
  pp2 <- posterior_predict(cs$ds, resp = "y2", newdata = nd2, ndraws = 20)
  ep2 <- posterior_epred(cs$ds, resp = "y2", newdata = nd2, ndraws = 20)
  expect_gt(mean(pp2[, 2]) / mean(ep2[, 2]), 0.9)
  expect_lt(mean(pp2[, 2]) / mean(ep2[, 2]), 1.1)
})
