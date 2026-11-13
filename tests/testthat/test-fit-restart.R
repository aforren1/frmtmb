# Fits that ended at a lower maximum with nothing to mark it
# (R/fit-restart.R): a prior far from the data's value of a coefficient,
# and a nonlinear body on a ridge that runs to a limit. Each defect test
# below was seen to fail on 0.81.0 (dev/nlstop-findings.md).

# A fit's messages, recorded and muffled, and the warnings `allowed`
# names, recorded and muffled by allow_warnings(); any other warning
# reaches testthat
fit_conds <- function(expr, allowed = character()) {
  got <- list()
  val <- allow_warnings(withCallingHandlers(expr,
    warning = function(w) got[[length(got) + 1L]] <<- w,
    message = function(m) {
      got[[length(got) + 1L]] <<- m
      invokeRestart("muffleMessage")
    }), allowed = allowed)
  list(value = val, conds = got)
}
has_class <- function(conds, cls) {
  any(vapply(conds, function(cnd) inherits(cnd, cls), NA))
}

# dev/review-nlstart-prior.R, rows apbx 1e6 1000: the response near 1e6,
# its slope 5e4, a normal(0, 1000) prior on the slope
conflict_data <- function(seed) {
  set.seed(seed)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$y <- 0.5 * d$x + rnorm(n)
  d$yl <- 1e6 * (1 + 0.05 * d$x) + rnorm(n, 0, 1e3)
  d
}
conflict_prior <- function(psd) {
  c(set_prior(sprintf("normal(0, %g)", psd), class = "b"),
    set_prior("normal(1e6, 1e6)", class = "Intercept"))
}

test_that("a prior far from the data's slope reaches the higher mode", {
  skip_on_cran()
  for (seed in c(6201, 6203)) {
    d <- conflict_data(seed)
    got <- fit_conds(frm(yl ~ x, data = d, prior = conflict_prior(1000)))
    fit <- got$value
    # the body spelling starts at the prior's location (prior_nl_starts())
    body <- fit_conds(frm(bf(yl ~ a + b * x, a + b ~ 1, nl = TRUE),
                          data = d,
                          prior = c(set_prior("normal(0, 1000)", nlpar = "b"),
                                    set_prior("normal(1e6, 1e6)",
                                              nlpar = "a"))))$value
    ll <- as.numeric(logLik_penalized(fit))
    ll_body <- as.numeric(logLik_penalized(body))
    # 0.81.0 stopped at the data's mode, about 470 below, at code 0
    expect_gte(ll, ll_body - 1e-8 * abs(ll_body))
    expect_lt(abs(fixef(fit)["x", "Estimate"]), 1e4)
    expect_true(has_class(got$conds, "frmtmb_prior_multistart"))
  }
})

test_that("a prior on the centered intercept far from the data moves it", {
  skip_on_cran()
  # the review's window: normal(0, 1e3) and 3e3 reach the prior's mode
  # from the default start already, from 2e4 the data's mode is the higher;
  # 1e4 is the case 0.81.0 stopped 3613.6 below the mode, at code 0
  d <- conflict_data(6201)
  pr <- set_prior("normal(0, 1e4)", class = "Intercept")
  got <- fit_conds(frm(yl ~ x, data = d, prior = pr))
  ref <- fit_conds(frm(yl ~ x, data = d, prior = pr,
                       start = list(beta = c(0, 5e4),
                                    betad = log(1e6))))$value
  ll <- as.numeric(logLik_penalized(got$value))
  ll_ref <- as.numeric(logLik_penalized(ref))
  expect_gte(ll, ll_ref - 1e-8 * abs(ll_ref))
  expect_true(has_class(got$conds, "frmtmb_prior_multistart"))
})

test_that("a prior the estimate agrees with takes no restart", {
  d <- conflict_data(6201)
  # the slope 5e4 is 0.05 prior scales from the prior's location
  got <- fit_conds(frm(yl ~ x, data = d, prior = conflict_prior(1e6)))
  expect_false(has_class(got$conds, "frmtmb_prior_multistart"))
  expect_gt(fixef(got$value)["x", "Estimate"], 4e4)
  # a given start is kept as given
  st <- list(beta = c(1e6, 5e4))
  got2 <- fit_conds(frm(yl ~ x, data = d, prior = conflict_prior(1000),
                        start = st))
  expect_false(has_class(got2$conds, "frmtmb_prior_multistart"))
  expect_gt(fixef(got2$value)["x", "Estimate"], 4e4)
})

test_that("the restart runs on a cold start with a free scale only", {
  d <- conflict_data(6201)
  fit <- fit_conds(frm(yl ~ x, data = d, prior = conflict_prior(1000)))$value
  calls <- 0L
  orig <- prior_conflict_escape
  local_mocked_bindings(prior_conflict_escape = function(...) {
    calls <<- calls + 1L
    orig(...)
  })
  # a refit is handed the fit's mode, as influence() and the bootstraps are
  rf <- fit_conds(refit(fit, d$yl * (1 + 1e-12)))
  expect_identical(calls, 0L)
  expect_false(has_class(rf$conds, "frmtmb_prior_multistart"))
  # a bernoulli likelihood has no dispersion to inflate: one mode, and the
  # restart is not run (a logistic GLMM paid twice its fit for it)
  set.seed(5)
  db <- data.frame(x = rnorm(300))
  db$y <- rbinom(300, 1, plogis(2 * db$x))
  fb <- fit_conds(frm(y ~ x, family = bernoulli(), data = db,
                      prior = set_prior("normal(0, 0.1)", class = "b")))
  expect_identical(calls, 1L)
  expect_false(has_class(fb$conds, "frmtmb_prior_multistart"))
  f <- fb$value
  # the slope is shrunk but still far from 0 in prior scales
  expect_gt(abs(fixef(f)["x", "Estimate"]) / 0.1, prior_conflict_z)
  local_mocked_bindings(optimize_obj = function(...) stop("restarted"))
  expect_identical(orig(f$obj, f$opt, f$frame, f$resolved_prior, f$control,
                        fit_outer_box(f), f$par_units), f$opt)
})

test_that("a zero-inflation probability is not a dispersion: no restart", {
  set.seed(5)
  d <- data.frame(x = rnorm(400))
  d$y <- ifelse(runif(400) < 0.3, 0L, rpois(400, exp(0.5 + d$x)))
  pr <- set_prior("normal(0, 0.05)", class = "b")
  got <- fit_conds(frm(y ~ x, family = zero_inflated_poisson(), data = d,
                       prior = pr))
  f <- got$value
  expect_false(has_class(got$conds, "frmtmb_prior_multistart"))
  # the slope is in conflict, and zi is the only distributional parameter
  expect_gt(abs(fixef(f)["x", "Estimate"]) / 0.05, prior_conflict_z)
  runs <- 0L
  local_mocked_bindings(optimize_obj = function(...) {
    runs <<- runs + 1L
    stop("restarted")
  })
  expect_identical(prior_conflict_escape(f$obj, f$opt, f$frame,
                                         f$resolved_prior, f$control,
                                         fit_outer_box(f), f$par_units),
                   f$opt)
  expect_identical(runs, 0L)
})

test_that("a conflicting prior under a one-mode likelihood keeps the fit", {
  set.seed(5)
  d <- data.frame(x = rnorm(300))
  d$y <- rpois(300, exp(0.5 + 1 * d$x))
  # the slope near 1 is 100 prior scales from 0; the posterior of a
  # log-concave likelihood under a normal prior has one mode
  got <- fit_conds(frm(y ~ x, family = poisson(), data = d,
                       prior = set_prior("normal(0, 0.01)", class = "b")))
  expect_false(has_class(got$conds, "frmtmb_prior_multistart"))
  ref <- fit_conds(frm(y ~ x, family = poisson(), data = d,
                       prior = set_prior("normal(0, 0.01)", class = "b"),
                       start = list(beta = unname(got$value$opt$par))))$value
  expect_equal(as.numeric(logLik_penalized(got$value)),
               as.numeric(logLik_penalized(ref)), tolerance = 1e-8)
})

# the review's fresh rows (dev/review-nlstart-fresh.R, body mm): every
# draw in the order it made them
mm_data <- function(seed, family) {
  set.seed(seed * 1000 + 1)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                  xq = round(runif(n, 0, 5), 3),
                  xb = round(runif(n, 0, 3000), 1),
                  g = factor(sample(1:15, n, TRUE)))
  d$ug <- rnorm(15, 0, 0.4)[d$g]
  lat <- 3 * d$xp / (1 + d$xp)
  d$y <- switch(family,
                bernoulli = rbinom(n, 1, plogis(lat - 2)),
                gaussian = lat + rnorm(n, 0, 0.5))
  d
}
mm_body <- bf(y ~ Vm * xp / (K + xp), Vm + K ~ 1, nl = TRUE)
# the two warnings a stop on the ridge can give
ridge_warn <- c("The data do not determine the coefficients",
                "Standard errors are not available")

test_that("a body on a ridge to its linear limit says so", {
  skip_on_cran()
  d <- mm_data(55, "bernoulli")
  got <- fit_conds(frm(mm_body, family = bernoulli(), data = d),
                   allowed = ridge_warn)
  fit <- got$value
  # 0.81.0 returned Vm 1.3e4 and K 1.8e5 at code 0 with no warning. The
  # fit ends at one of two stops by rounding: at the linear limit, where
  # the limit warning is given, or 0.0065 below it, where the
  # standard-error check's own warning is (17 and 3 of 20 ulp draws,
  # dev/nlstop-limit.R; OpenBLAS 0.3.26 takes the second). Either warns,
  # and the check itself finds the ridge from either
  expect_true(any(vapply(got$conds, inherits, NA, "warning")))
  expect_false(is.null(nl_limit_message(fit)))
  if (has_class(got$conds, "frmtmb_nl_limit")) {
    lin <- frm(y ~ 0 + xp, family = bernoulli(), data = d)
    ll_lin <- as.numeric(logLik(lin))
    # at the limit: the linear model's value to grad_tol
    expect_lt(abs(as.numeric(logLik(fit)) - ll_lin), fit$control$grad_tol)
  }
  # and so on every rounding of the covariate a few ulps apart
  for (k in 1:4) {
    dk <- d
    dk$xp <- dk$xp * (1 + 2e-15 * k)
    gk <- fit_conds(frm(mm_body, family = bernoulli(), data = dk),
                    allowed = ridge_warn)
    expect_true(any(vapply(gk$conds, inherits, NA, "warning")))
  }
})

test_that("an identified body, or a coefficient near 0, gives no limit", {
  # no warning is allowed: any would reach testthat
  d <- mm_data(55, "gaussian")
  got <- fit_conds(frm(mm_body, data = d))
  expect_false(has_class(got$conds, "frmtmb_nl_limit"))
  set.seed(3)
  d2 <- data.frame(x = rnorm(300))
  d2$y <- 2 + rnorm(300, 0, 0.5)
  got2 <- fit_conds(frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
                        data = d2))
  expect_false(has_class(got2$conds, "frmtmb_nl_limit"))
})
