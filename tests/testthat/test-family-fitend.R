# What a family says where its fit lands, and the family gradients it
# needs there: the GEV's non-regular range and its far lower tail, the
# Frechet at nu = 1, the cs() ordinal-mixture wall, a negbinomial at
# its poisson limit in either spelling, and a start where a prior's
# density is 0 (lane famend, dev/famend-findings.md).

fe_conds <- function(expr) {
  w <- list()
  m <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w[[length(w) + 1L]] <<- x
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m[[length(m) + 1L]] <<- x
    invokeRestart("muffleMessage")
  })
  list(value = v, warnings = w, messages = m)
}
fe_has <- function(conds, cls) {
  any(vapply(conds, function(x) inherits(x, cls), NA))
}
# every warning fe_conds() caught is one the test names, by class or by
# a fixed piece of its text, and any other fails the test
fe_only <- function(r, allowed) {
  ok <- vapply(r$warnings, function(x) {
    any(vapply(allowed, function(a) {
      inherits(x, a) || grepl(a, conditionMessage(x), fixed = TRUE)
    }, NA))
  }, NA)
  expect_true(all(ok), info = paste(vapply(r$warnings[!ok], conditionMessage,
                                           ""), collapse = " || "))
}

fe_rgev <- function(n, mu, s, xi) {
  mu + s * ((-log(stats::runif(n)))^(-xi) - 1) / xi
}

test_that("a GEV fit with xi below -0.5 says its maximum is not regular", {
  # PIN: 0.82.0 returned the standard errors of these fits in silence
  # (dev/famend-repro.R, case gev: xi -0.663 and -0.936 on seed 4)
  for (xi0 in c(-0.63, -0.89)) {
    set.seed(4)
    d <- data.frame(y = fe_rgev(300, 0, 1, xi0))
    r <- fe_conds(frm(y ~ 1, family = gen_extreme_value(), data = d))
    expect_length(r$warnings, 0L)
    expect_true(fe_has(r$messages, "frmtmb_gev_xi_nonregular"),
                info = paste("xi0", xi0))
    # the profile call the message gives runs as written (it named no
    # parm at first, and confint() needs one for a profile)
    msg <- conditionMessage(r$messages[[1L]])
    parm <- sub('.*parm = "([^"]+)".*', "\\1", msg)
    expect_identical(parm, "xi_(Intercept)")
    ci <- allow_warnings(confint(r$value, parm = parm, method = "profile"),
                         "collapsing to unique 'x' values")
    expect_true(is.finite(ci[1L, "upr"]))
    xi <- frm_linpred(r$value, dpar = "xi", type = "response")[1]
    expect_lt(xi, -0.5)
    expect_gt(xi, -1 + 1e-3)
  }
  # the guard absent: xi -0.45 ends at -0.473, inside the regular range
  set.seed(4)
  d <- data.frame(y = fe_rgev(300, 0, 1, -0.45))
  r <- fe_conds(frm(y ~ 1, family = gen_extreme_value(), data = d))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 0L)
  expect_gt(frm_linpred(r$value, dpar = "xi", type = "response")[1], -0.5)
})

test_that("a GEV row far below its support has a finite gradient", {
  # PIN: on 0.82.0 an outside row 19.47 to 20.26 scale units below the
  # location had a NaN gradient at xi from 0.06 to 10: exp(-w), up to
  # exp(700), times the 1e15 slope of log1p_ratio() at the clamp
  # overflowed before the clamp's zero met it (dev/famend-gevbit.R: 16
  # of 1200 sampled rows at each xi). At xi 0.44 to 0.84 the clamped z
  # rounded 7 ulps above the clamp, past the ratio's hold, and 490 of
  # 901 rows stayed NaN after a first fix
  # (dev/reviews/2026-11-13-famend-short.md)
  fam <- gen_extreme_value()
  xs <- seq(-19.4, -20.3, by = -0.01)
  for (xi in c(0.06, 0.5, 0.7, 1, 10)) {
    p <- c(0, 0, log1p(xi))
    bad <- 0L
    nan_h <- 0L
    for (y in xs) {
      tp <- RTMB::MakeTape(function(p) {
        fam$lpdf(y, list(mu = p[1], sigma = exp(p[2]), xi = expm1(p[3]),
                         .eta_sigma = p[2], .eta_xi = p[3]), list())
      }, p)
      bad <- bad + any(!is.finite(tp$jacobian(p)))
      # a Hessian entry near h = -700 is about exp(700) 34.5^2 x^2, past
      # the doubles' end, so it may be Inf; it must not be NaN
      nan_h <- nan_h + any(is.nan(tp$jacfun()$jacobian(p)))
    }
    expect_identical(bad, 0L, info = paste("xi", xi))
    expect_identical(nan_h, 0L, info = paste("xi", xi))
  }
})

fe_frechet_data <- function(seed, n = 300, nu0 = 1.1) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / nu0) * stats::rexp(n)^(-1 / nu0)
  d
}

test_that("a Frechet fit that runs to nu = 1 names that edge", {
  # PIN: on 0.82.0 seed 9031 ran to nu = 1 + 4.9e-9 with the mean
  # intercept at 17.38, and only the standard-error warning was given
  # (dev/famend-repro.R, case frechet)
  d31 <- fe_frechet_data(9031)
  r <- fe_conds(frm(y ~ x, family = frechet(), data = d31))
  expect_true(fe_has(r$warnings, "frmtmb_frechet_nu_edge"))
  # the SE warning it explains is not repeated, and nothing else is said
  fe_only(r, "frmtmb_frechet_nu_edge")
  f <- r$value
  nu <- frm_linpred(f, dpar = "nu", type = "response")[1]
  expect_lt(nu - 1, 1e-3)
  # PIN: the check's probe toward the edge is better than the stop, and
  # TMB kept it as its best point, which sdreport() and the profile read
  # (dev/reviews/2026-11-13-famend.md, F2). The fit's state is the stop,
  # and every SE and profile is the one the fit gives without the check
  env <- f$obj$env
  expect_identical(unname(env$last.par.best), unname(f$opt$par))
  expect_identical(env$value.best, f$opt$objective)
  local({
    local_mocked_bindings(frechet_fit_check = function(fit, resp) {
      invisible(NULL)
    }, .package = "frmtmb")
    f_off <- allow_warnings(frm(y ~ x, family = frechet(), data = d31),
                            "Standard errors are not available",
                            require = "Standard errors are not available")
    expect_identical(vcov(f, full = TRUE), vcov(f_off, full = TRUE))
    expect_identical(confint(f, parm = "x", method = "profile"),
                     confint(f_off, parm = "x", method = "profile"))
  })
  # the guard absent: seed 9001 has a maximum at nu - 1 = 0.116, and
  # moving toward the edge with the scale held loses likelihood there
  r0 <- fe_conds(frm(y ~ x, family = frechet(), data = fe_frechet_data(9001)))
  fe_only(r0, character(0))
  expect_length(r0$messages, 0L)
  f0 <- r0$value
  lp <- f0$frame[["linpreds"]][[frmtmb:::linpred_key(
    f0$spec$responses[[1]]$resp_name, "nu")]]
  expect_lt(frmtmb:::frechet_ridge_gain(f0, f0$spec$responses[[1]]$resp_name,
                                        lp), 0)
})

test_that("a cs() ordinal mixture at its crossing wall names the wall", {
  skip_on_cran()
  mk <- function(seed, n = 300) {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
    cls <- stats::rbinom(n, 1, 0.4)
    lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
    d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
    d
  }
  fit_cs <- function(d) {
    fe_conds(frm(bf(y ~ x + cs(z)), family = mixture(cumulative(), sratio()),
                 data = d))
  }
  wall <- function(r) {
    any(grepl("cs() effect has closed a category",
              vapply(r$warnings, conditionMessage, ""), fixed = TRUE))
  }
  # brms's own warning on cs() for these families
  brms_cs <- "Category specific effects for this family"
  # PIN: on 0.82.0 seed 3 stopped at false convergence with a row's
  # category 1.4e-14 wide, and no warning named it (dev/famend-repro.R,
  # case csmix)
  r <- fit_cs(mk(3))
  expect_true(wall(r))
  fe_only(r, c(brms_cs, "frmtmb_conv_code", "frmtmb_conv_gradient",
               "cs() effect has closed a category"))
  # the guard absent: seed 4 converges, its narrowest category 0.509;
  # its component 2 is a step function, which the same warning names
  r4 <- fit_cs(mk(4))
  expect_false(wall(r4))
  fe_only(r4, c(brms_cs, "a step function of its predictors"))
})

test_that("negbinomial and negbinomial2 both name the poisson limit", {
  # both spellings stop on the ridge where the likelihood is flat, on its
  # boundary (dev/famend-repro.R, case nbmirror; backlog item of lane
  # nb2dm). Where the stop fell on that ridge was rounding, through the
  # noise of RTMB::dnbinom_robust()'s derivative in the shape, and 22 of
  # 84 roundings stopped at code 1; with the exact derivative all 84
  # stop at code 0 with the message (dev/nbrobust-findings.md). The
  # log-likelihood is the limit's wherever the stop is: the stops of
  # seeds 4 and 6 differ by 2.7e-6 at most, below 1e-8 of its 555
  set.seed(4)
  d <- data.frame(x = stats::rnorm(300))
  d$y <- stats::rpois(300, exp(1 + 0.3 * d$x))
  r1 <- fe_conds(frm(y ~ x, family = negbinomial(), data = d))
  r2 <- fe_conds(frm(y ~ x, family = brmsfamily("negbinomial2"), data = d))
  for (r in list(r1, r2)) {
    expect_true(fe_has(r$messages, "frmtmb_limit_fit") ||
                  fe_has(r$messages, "frmtmb_boundary_fit"))
    fe_only(r, character(0))
  }
  l1 <- as.numeric(logLik(r1$value))
  l2 <- as.numeric(logLik(r2$value))
  expect_lt(abs(l1 - l2) / abs(l1), 1e-8)
})

test_that("a negbinomial fit started far along its poisson limit stops there", {
  # PIN: started at a log shape of 30, 15 of 21 roundings of these data
  # warned "Large maximum absolute gradient" (0.191 at the shape) and
  # "false convergence" on 0.83.0, and all 21 under check_se = "ignore".
  # That gradient was the noise of RTMB::dnbinom_robust()'s derivative;
  # the true one is about 1e-13, so the start is already the limit
  # (dev/nbrobust-start.R)
  set.seed(2)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.4 * d$x))
  r <- fe_conds(frm(y ~ x, family = negbinomial(), data = d,
                    start = list(betad = 30)))
  expect_true(fe_has(r$messages, "frmtmb_limit_fit") ||
                fe_has(r$messages, "frmtmb_boundary_fit"))
  fe_only(r, character(0))
  r <- fe_conds(frm(y ~ x, family = negbinomial(), data = d,
                    start = list(betad = 30),
                    control = frmtmb_control(check_se = "ignore")))
  fe_only(r, character(0))
})

test_that("a start where a prior's density is 0 moves inside its support", {
  # PIN: on 0.82.0 gamma(2, 1) on the coefficients of shape ~ 0 + fac,
  # which start at 0, left the fit at its start with logLik NaN, with or
  # without lb = 0 (dev/famend-probe-gprior.R)
  set.seed(20261021)
  n <- 80
  d <- data.frame(x = stats::rnorm(n), fac = factor(rep(c("a", "b"), 40)))
  d$yg <- stats::rgamma(n, shape = 3, rate = 3 / exp(0.5 + 0.4 * d$x))
  fo <- bf(yg ~ x, shape ~ 0 + fac)
  pr <- set_prior("gamma(2, 1)", class = "b", dpar = "shape")
  st <- allow_warnings(par_template(fo, data = d, family = Gamma("log"),
                                    prior = pr),
                       allowed = "lower bounded prior")
  expect_identical(unname(st$betad), c(2, 2))
  f1 <- allow_warnings(frm(fo, family = Gamma("log"), data = d, prior = pr),
                       allowed = "lower bounded prior")
  f2 <- frm(fo, family = Gamma("log"), data = d,
            prior = set_prior("gamma(2, 1)", class = "b", dpar = "shape",
                              lb = 0))
  expect_identical(f1$opt$convergence, 0L)
  expect_true(is.finite(as.numeric(logLik(f1))))
  expect_true(all(f1$estimates$betad > 0))
  expect_lt(abs(f1$opt$objective - f2$opt$objective) / abs(f2$opt$objective),
            1e-8)
  # the guard absent: a prior whose density is positive at the start
  # leaves the start alone
  st0 <- allow_warnings(
    par_template(fo, data = d, family = Gamma("log"),
                 prior = set_prior("exponential(1)", class = "b",
                                   dpar = "shape")),
    allowed = "lower bounded prior")
  expect_identical(unname(st0$betad), c(0, 0))
})
