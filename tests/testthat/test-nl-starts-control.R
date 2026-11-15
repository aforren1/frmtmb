# frmtmb_control(nl_starts = ): FALSE turns off two automatic restarts of
# a nonlinear fit, the sign probe and the restart beyond a pole of the
# body; TRUE (the default) keeps them (user decision of 2026-11-13). The
# restart from a saddle point is outside the switch.

nls_conds <- function(expr, allowed = character()) {
  got <- list()
  val <- allow_warnings(withCallingHandlers(expr,
    warning = function(w) got[[length(got) + 1L]] <<- w,
    message = function(m) {
      got[[length(got) + 1L]] <<- m
      invokeRestart("muffleMessage")
    }), allowed = allowed)
  list(value = val, conds = got)
}
nls_has <- function(conds, cls) {
  any(vapply(conds, function(cnd) inherits(cnd, cls), NA))
}
# the classes of a restart's message or warning
nls_restart <- c("frmtmb_nl_multistart", "frmtmb_nl_multistart_unsettled",
                 "frmtmb_nl_pole_start", "frmtmb_nl_pole")
# the warnings a stop on a ridge or between poles can give
nls_ridge <- c("The data do not determine the coefficients",
               "Standard errors are not available")
nls_mm <- function(seed, family) {
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
                binomial = rbinom(n, 10, plogis(lat - 2)))
  d
}

test_that("nl_starts is a flag, TRUE by default", {
  expect_true(frmtmb_control()$nl_starts)
  expect_false(frmtmb_control(nl_starts = FALSE)$nl_starts)
  expect_error(frmtmb_control(nl_starts = "no"), "nl_starts")
  expect_error(frmtmb_control(nl_starts = c(TRUE, FALSE)), "nl_starts")
})

test_that("nl_starts = FALSE keeps the fit between the poles; TRUE leaves", {
  skip_on_cran()
  # binomial Vm * xp / (K + xp), xp in 0.1 to 5, seeds 11 and 44: the fit
  # from the default start stops between the poles of rows inside the
  # data, and the pole restart takes it beyond them, to a higher maximum
  fm <- bf(y | trials(nt) ~ Vm * xp / (K + xp), Vm + K ~ 1, nl = TRUE)
  for (seed in c(11, 44)) {
    d <- nls_mm(seed, "binomial")
    on <- nls_conds(frm(fm, family = binomial(), data = d))
    off <- nls_conds(frm(fm, family = binomial(), data = d,
                         control = frmtmb_control(nl_starts = FALSE)),
                     allowed = nls_ridge)
    tol <- on$value$control$grad_tol
    expect_true(nls_has(on$conds, "frmtmb_nl_pole_start"))
    expect_lt(on$value$estimates$beta[[2L]], -max(d$xp))
    expect_false(nls_has(off$conds, nls_restart))
    expect_gt(off$value$estimates$beta[[2L]], -max(d$xp))
    expect_lt(as.numeric(logLik(off$value)),
              as.numeric(logLik(on$value)) - tol)
  }
})

test_that("nl_starts = FALSE keeps a given start's fit without the check", {
  skip_on_cran()
  # bernoulli seed 55 from K -2.5 stops at a pole inside the data; with
  # the restarts on, the run beyond the poles says that this is not the
  # maximum (frmtmb_nl_pole); with them off nothing is run or said
  d <- nls_mm(55, "bernoulli")
  body <- bf(y ~ Vm * xp / (K + xp), Vm + K ~ 1, nl = TRUE)
  st <- list(beta = c(-2, -2.5))
  on <- nls_conds(frm(body, family = bernoulli(), data = d, start = st),
                  allowed = c("A start beyond the poles", nls_ridge))
  off <- nls_conds(frm(body, family = bernoulli(), data = d, start = st,
                       control = frmtmb_control(nl_starts = FALSE)),
                   allowed = nls_ridge)
  expect_true(nls_has(on$conds, "frmtmb_nl_pole"))
  expect_false(nls_has(off$conds, nls_restart))
  # the given start's fit is kept either way
  expect_identical(on$value$opt$par, off$value$opt$par)
})

test_that("nl_starts = FALSE switches off the sign probe", {
  # b1 * (1 - exp(-b2 * x)) with its maximum at both signs negative: the
  # probe finds it from the default start, and without the probe the fit
  # stops short, as under options(frmtmb.nl_sign_probe = FALSE)
  set.seed(101)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  d$y <- 1.5 * exp(0.6 * d$x) + 0.3 * d$z + rnorm(n, 0, 0.5)
  body <- bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE)
  # without the probe the fit stops on the ridge to the linear limit, where
  # nlminb reports false convergence
  stop8 <- c(nls_ridge, "on a ridge", "Optimizer did not report convergence")
  on <- nls_conds(frm(body, data = d), allowed = nls_ridge)
  off <- nls_conds(frm(body, data = d,
                       control = frmtmb_control(nl_starts = FALSE)),
                   allowed = stop8)
  opt <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
    nls_conds(frm(body, data = d), allowed = stop8))
  expect_true(nls_has(on$conds, "frmtmb_nl_multistart"))
  expect_false(nls_has(off$conds, nls_restart))
  ll_on <- as.numeric(logLik(on$value))
  expect_gt(ll_on - as.numeric(logLik(off$value)), 1e-3 * abs(ll_on))
  expect_identical(off$value$opt$par, opt$value$opt$par)
})
