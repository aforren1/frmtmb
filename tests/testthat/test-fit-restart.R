# Fits that ended at a lower maximum with nothing to mark it
# (R/fit-restart.R): a prior far from the data's value of a coefficient,
# and a nonlinear body on a ridge that runs to a limit. Each defect test
# below was seen to fail on 0.81.0 (dev/nlstop-findings.md); those from a
# bending ridge on were seen to fail on 0.82.0 (dev/nlstop2-findings.md).
# Then a body's poles among the rows, and in a hole of the data (seen to
# fail on 0.83.0, dev/nlband-findings.md).

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
  d$nt <- 10L
  d$y <- switch(family,
                bernoulli = rbinom(n, 1, plogis(lat - 2)),
                binomial = rbinom(n, 10, plogis(lat - 2)),
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

test_that("a ridge that bends is read after a step across it", {
  # a * (exp(b * x) - 1) on linear data: 0.82.0 stopped at a 1.40e4, b
  # 1.44e-4, code 0, with no warning, 0.0033 below the ridge's own value at
  # (1000 a, b / 1000) (review of lane nlstop, MINOR 1)
  set.seed(2)
  d <- data.frame(x = runif(1000, 0, 3))
  d$y <- 2 * d$x + rnorm(1000, 0, 0.5)
  # by rounding the fit ends on the ridge at code 0 with the limit
  # warning, at code 0 with the gradient warning, or further along at code
  # 1 with the optimizer's warnings (dev/nlstop2-ulp.R explin; review of
  # lane nlstop2, F1); every ending warns, and the check itself finds the
  # ridge from each
  got <- fit_conds(frm(bf(y ~ a * (exp(b * x) - 1), a + b ~ 1, nl = TRUE),
                       data = d, start = list(beta = c(2, 0.5))),
                   allowed = c("The data do not determine the coefficients",
                               "Optimizer did not report convergence",
                               "Large maximum absolute gradient"))
  expect_true(any(vapply(got$conds, inherits, NA, "warning")))
  expect_false(is.null(nl_limit_message(got$value)))
})

test_that("a fit between the poles of its body is taken beyond them", {
  skip_on_cran()
  # binomial Vm * xp / (K + xp), xp in 0.1 to 5: 0.82.0 stopped at code 0
  # with no warning at K -0.1825 (seed 11), between the poles of rows
  # inside the data, and at K -0.1103 (seed 44), beside the pole of the
  # smallest xp; a start at K -6 reaches the maxima at K -6.0152 and
  # -6.4531, 4.32 and 0.92 higher (review of lane nlstop, MINOR 3)
  for (seed in c(11, 44)) {
    d <- mm_data(seed, "binomial")
    fm <- bf(y | trials(nt) ~ Vm * xp / (K + xp), Vm + K ~ 1, nl = TRUE)
    ref <- frm(fm, family = binomial(), data = d,
               start = list(beta = c(-0.1, -6)))
    tol <- ref$control$grad_tol
    got <- fit_conds(frm(fm, family = binomial(), data = d))
    fit <- got$value
    expect_true(has_class(got$conds, "frmtmb_nl_pole_start"))
    expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), tol)
    expect_lt(fit$estimates$beta[[2L]], -max(d$xp))
    # every rounding of the covariate a few ulps apart
    for (k in 1:4) {
      dk <- d
      dk$xp <- dk$xp * (1 + 2e-15 * k)
      gk <- fit_conds(frm(fm, family = binomial(), data = dk))
      expect_lt(as.numeric(logLik(ref)) - as.numeric(logLik(gk$value)), tol)
    }
  }
})

test_that("a given start between the poles keeps its fit, with a warning", {
  skip_on_cran()
  # bernoulli seed 55 from K -2.5: 0.82.0 stopped at code 0 with no
  # warning at K -1.00, a pole inside the data, 1.38 below the linear
  # limit the start beyond the poles reaches
  d <- mm_data(55, "bernoulli")
  got <- fit_conds(frm(mm_body, family = bernoulli(), data = d,
                       start = list(beta = c(-2, -2.5))),
                   allowed = "A start beyond the poles")
  fit <- got$value
  expect_true(has_class(got$conds, "frmtmb_nl_pole"))
  K <- fit$estimates$beta[[2L]]
  expect_true(K > -max(d$xp) && K < -min(d$xp))
})

test_that("a body clear of its poles takes no start beyond them", {
  # gaussian seed 55: K 0.92, 1.05 from the nearest pole, whose gaps are
  # 0.012; no warning is allowed, so any would reach testthat
  d <- mm_data(55, "gaussian")
  got <- fit_conds(frm(mm_body, data = d))
  fit <- got$value
  expect_false(has_class(got$conds, "frmtmb_nl_pole_start"))
  lp <- Filter(function(q) !is.null(q$nl_body), fit$frame$linpreds)[[1L]]
  wl <- nl_pole_walls(fit$frame, lp, fit$estimates, "beta", 2L,
                      nl_denoms(lp$nl_body))
  # the walls of K are read, at -xp, and K is clear of them
  expect_length(wl, 1L)
  expect_equal(c(wl[[1L]]$lo, wl[[1L]]$hi), -rev(range(d$xp)),
               tolerance = 1e-8)
  expect_false(nl_pole_pressed(wl[[1L]]))
  expect_null(nl_pole_escape(function(...) stop("no refit"), fit$frame, fit))
  # walls on both sides of the fit, 1 + b * x with x of both signs: the
  # denominator is positive on every row, and the fit is 10.4 median gaps
  # from the nearest wall (dev/nlstop2-gaps.R, seed 1)
  set.seed(1)
  d2 <- data.frame(x = rnorm(300))
  d2$y <- 2 / (1 + 0.1 * d2$x) + rnorm(300, 0, 0.3)
  got2 <- fit_conds(frm(bf(y ~ a / (1 + b * x), a + b ~ 1, nl = TRUE),
                        data = d2))
  f2 <- got2$value
  lp2 <- Filter(function(q) !is.null(q$nl_body), f2$frame$linpreds)[[1L]]
  wl2 <- nl_pole_walls(f2$frame, lp2, f2$estimates, "beta", 2L,
                       nl_denoms(lp2$nl_body))
  expect_true(wl2[[1L]]$lo < f2$estimates$beta[[2L]] &&
                wl2[[1L]]$hi > f2$estimates$beta[[2L]])
  expect_false(nl_pole_pressed(wl2[[1L]]))
  expect_null(nl_pole_escape(function(...) stop("no refit"), f2$frame, f2))
})

# y = 1 / (x - 2) + noise with no rows within 0.3 of 2: the pole of a / (x
# + b) belongs in that gap of the data, at b -2 (review of lane nlstop2,
# F4); `inv` = 1 / (x + inv) is the truth beyond the data instead
band_data <- function(seed, inv = -2) {
  set.seed(seed)
  x <- runif(400, 0, 4)
  x <- x[abs(x - 2) > 0.3]
  data.frame(x = x, y = 1 / (x + inv) + rnorm(length(x), 0, 0.2))
}
band_body <- bf(y ~ a / (x + b), a + b ~ 1, nl = TRUE)
# the walls of b a fit reads, and the hole of the data between them
band_walls <- function(fit) {
  lp <- Filter(function(q) !is.null(q$nl_body), fit$frame$linpreds)[[1L]]
  nl_pole_walls(fit$frame, lp, fit$estimates, "beta", 2L,
                nl_denoms(lp$nl_body))
}

test_that("a fit clear of its poles is started in a hole of the data", {
  # 0.83.0 stopped at b 0.39, 0.48, 0.37 with every denominator positive,
  # code 0 and no warning, 592, 592 and 646 below the fit from b -1.9
  for (seed in 1:3) {
    d <- band_data(seed)
    ref <- frm(band_body, data = d, start = list(beta = c(1, -1.9)))
    tol <- ref$control$grad_tol
    got <- fit_conds(frm(band_body, data = d))
    fit <- got$value
    expect_true(has_class(got$conds, "frmtmb_nl_pole_start"))
    expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))), tol)
    b <- fit$estimates$beta[[2L]]
    expect_true(b > -min(d$x[d$x > 2]) && b < -max(d$x[d$x < 2]))
  }
  # every rounding of the response a few ulps apart
  d <- band_data(1)
  ref <- frm(band_body, data = d, start = list(beta = c(1, -1.9)))
  for (k in 1:4) {
    dk <- d
    dk$y <- dk$y * (1 + 2e-15 * k)
    gk <- fit_conds(frm(band_body, data = dk))
    expect_true(has_class(gk$conds, "frmtmb_nl_pole_start"))
    expect_lt(as.numeric(logLik(ref)) - as.numeric(logLik(gk$value)),
              ref$control$grad_tol)
  }
})

test_that("a given start clear of the poles warns of the hole", {
  d <- band_data(1)
  got <- fit_conds(frm(band_body, data = d,
                       start = list(beta = c(-0.5, 0.4))),
                   allowed = "a gap that the rows' poles leave")
  fit <- got$value
  expect_true(has_class(got$conds, "frmtmb_nl_pole"))
  # the caller's start is kept: the fit is still clear of the poles
  expect_gt(fit$estimates$beta[[2L]], -min(d$x))
})

test_that("a hole whose start is worse, or no hole, keeps the fit", {
  # no warning is allowed: any would reach testthat. The escape catches an
  # error of its refit, so a call is recorded, not only refused
  called <- new.env()
  refit <- function(...) {
    called$refit <- TRUE
    stop("the escape refit a losing start")
  }
  for (seed in 1:2) {
    # the truth beyond the data: the hole start runs and loses
    d <- band_data(seed, inv = 0.5)
    got <- fit_conds(frm(band_body, data = d))
    fit <- got$value
    expect_false(has_class(got$conds, "frmtmb_nl_pole_start"))
    expect_gt(fit$estimates$beta[[2L]], -min(d$x))
    wl <- band_walls(fit)
    expect_length(wl, 1L)
    expect_identical(nrow(wl[[1L]]$holes), 1L)
    expect_null(nl_pole_escape(refit, fit$frame, fit))
    expect_null(called$refit)
    # every row kept: no hole, so no start
    set.seed(seed)
    d2 <- data.frame(x = runif(400, 0, 4))
    d2$y <- 1 / (d2$x + 0.5) + rnorm(400, 0, 0.2)
    got2 <- fit_conds(frm(band_body, data = d2))
    expect_false(has_class(got2$conds, "frmtmb_nl_pole_start"))
    wl2 <- band_walls(got2$value)
    expect_length(wl2, 1L)
    expect_null(wl2[[1L]]$holes)
  }
  # nl_starts = FALSE leaves the band data's fit where it stopped
  d <- band_data(1)
  off <- fit_conds(frm(band_body, data = d,
                       control = frmtmb_control(nl_starts = FALSE)))
  expect_false(has_class(off$conds, "frmtmb_nl_pole_start"))
  expect_gt(off$value$estimates$beta[[2L]], -min(d$x))
  # and from that fit, with the switch on, the hole start wins, so the
  # recorder above can see a refit
  f0 <- off$value
  f0$control$nl_starts <- TRUE
  expect_null(nl_pole_escape(refit, f0$frame, f0))
  expect_true(isTRUE(called$refit))
})

test_that("a hole is a gap much wider than the gaps around it", {
  # two even runs of 50 walls: the gap between them is the one hole
  w <- c(seq(0, 1, length.out = 50), seq(2, 3, length.out = 50))
  h <- nl_pole_holes(w)
  expect_identical(nrow(h), 1L)
  expect_equal(unname(h[1L, ]), c(1, 2))
  # an even run has none, and neither has a wide gap with fewer than
  # nl_hole_side walls on one side of it
  expect_null(nl_pole_holes(seq(0, 3, length.out = 100)))
  expect_null(nl_pole_holes(c(seq(0, 0.08, length.out = 5),
                              seq(1, 3, length.out = 100))))
  expect_null(nl_pole_holes(c(seq(0, 1, length.out = 10),
                              seq(2, 3, length.out = 10))))
  # three holes: the nl_hole_max widest, widest first
  w3 <- c(seq(0, 1, length.out = 50), seq(2, 3, length.out = 50),
          seq(6, 7, length.out = 50), seq(9, 10, length.out = 50))
  h3 <- nl_pole_holes(w3)
  expect_identical(nrow(h3), nl_hole_max)
  expect_equal(unname(h3[, "lo"]), c(3, 7))
})

test_that("a hole run that rises without end is not kept", {
  skip_on_cran()
  # y = 1 for x > 2 on the band data: a pole in the hole separates every
  # row, and the likelihood rises without end as a grows. The first hole
  # build kept the run on seeds 2 to 10: a 1052 to 1237, logLik about
  # -1e-268, standard errors about 1e133, code 0, no warning, with the
  # message (review of lane nlband, item 1); 0.83.0 warned on all 10
  stop_warn <- c("Optimizer did not report convergence",
                 "Large maximum absolute gradient",
                 "Standard errors are not available")
  for (seed in 2:4) {
    d <- band_data(seed)
    d$y <- as.integer(d$x > 2)
    got <- fit_conds(frm(band_body, family = bernoulli(), data = d),
                     allowed = stop_warn)
    expect_false(has_class(got$conds, "frmtmb_nl_pole_start"))
    expect_true(any(vapply(got$conds, inherits, NA, "warning")))
    # the refused fit is higher, so frm() says why it is not returned
    expect_true(has_class(got$conds, "frmtmb_nl_pole_refused"))
  }
})

test_that("a hole run starts the rest of the model at the default start", {
  skip_on_cran()
  # student band data: the clear fit has log nu at its limit, 19.7, and
  # hole runs that kept it there ended 16.2, 20.3 and 44.0 below the
  # maximum, with the boundary message of the clear fit (review of lane
  # nlband, item 2)
  for (seed in 1:3) {
    set.seed(seed)
    x <- band_data(seed)$x
    d <- data.frame(x = x, y = 1 / (x - 2) + 0.2 * rt(length(x), 3))
    ref <- frm(band_body, data = d, family = student(),
               start = list(beta = c(1, -2)),
               control = frmtmb_control(nl_starts = FALSE))
    off <- suppressMessages(frm(band_body, data = d, family = student(),
                                control = frmtmb_control(nl_starts = FALSE)))
    got <- fit_conds(frm(band_body, data = d, family = student()))
    expect_true(has_class(got$conds, "frmtmb_nl_pole_start"))
    expect_false(has_class(got$conds, "frmtmb_boundary_fit"))
    # within a millionth of the gap the hole start closes (481 to 521)
    gap <- as.numeric(logLik(ref)) - as.numeric(logLik(off))
    expect_lt(abs(as.numeric(logLik(got$value)) - as.numeric(logLik(ref))),
              1e-6 * gap)
  }
})

test_that("a kept hole fit says its own boundary message", {
  skip_on_cran()
  # student() on the gaussian band data, seed 2: the maximum has nu at its
  # limit; the fit from (1, -2) says so, and the kept hole fit said
  # nothing, its refit's messages suppressed (short check of lane nlband,
  # item 1)
  d <- band_data(2)
  ref <- fit_conds(frm(band_body, data = d, family = student(),
                       start = list(beta = c(1, -2)),
                       control = frmtmb_control(nl_starts = FALSE)))
  got <- fit_conds(frm(band_body, data = d, family = student()))
  expect_true(has_class(got$conds, "frmtmb_nl_pole_start"))
  expect_true(has_class(ref$conds, "frmtmb_boundary_fit"))
  expect_true(has_class(got$conds, "frmtmb_boundary_fit"))
})

test_that("a hole fit is kept where the data place its pole loosely", {
  skip_on_cran()
  # n 60 with the narrowest hole about 2 the rule flags, y = 0.2 / (x - 2)
  # + N(0, 2), seed 13: the hole's refit is a maximum 1.46 above the first
  # fit with a standard error of b 0.29, above half the hole's width,
  # 0.27, and far below the span of the poles, 3.9; the round-1 rule
  # refused it in silence (short check of lane nlband, item 2)
  for (h in seq(0.2, 1.6, by = 0.02)) {
    set.seed(13)
    x <- runif(2000, 0, 4)
    x <- x[abs(x - 2) > h][seq_len(60)]
    hs <- nl_pole_holes(sort(unique(-x)))
    if (!is.null(hs) && any(hs[, "lo"] < -2 & hs[, "hi"] > -2)) break
  }
  set.seed(513)
  d <- data.frame(x = x, y = 0.2 / (x - 2) + rnorm(60, 0, 2))
  off <- suppressMessages(frm(band_body, data = d,
                              control = frmtmb_control(nl_starts = FALSE)))
  got <- fit_conds(frm(band_body, data = d))
  expect_true(has_class(got$conds, "frmtmb_nl_pole_start"))
  b <- got$value$estimates$beta[[2L]]
  expect_true(b > -min(d$x[d$x > 2]) && b < -max(d$x[d$x < 2]))
  expect_gt(as.numeric(logLik(got$value)), as.numeric(logLik(off)))
})
