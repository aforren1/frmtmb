# logLik() under a prior is the log-likelihood of the data at the MAP
# estimates; logLik_penalized() is the value the fit maximized, with the
# log prior. Through frmtmb 0.79.1 logLik() was the penalized value, so
# a normal(0, 100) prior on one slope lowered it by log(100 * sqrt(2 pi))
# = 5.52 at the same estimates, and AIC(m0, m1) compared the two with no
# word (dev/loglik-findings.md, round r080).

lp_data <- function(seed = 4101L, n = 120L, ng = 12L) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(rep(seq_len(ng), each = n / ng)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(ng, 0, 0.7)[d$g] + stats::rnorm(n)
  d
}

lp_prior <- function() {
  set_prior("normal(0, 0.2)", class = "b") +
    set_prior("student_t(3, 0, 1)", class = "sd") +
    set_prior("normal(2, 0.5)", class = "Intercept")
}

test_that("logLik() under a prior is the closed-form likelihood at the MAP", {
  set.seed(7)
  d <- data.frame(x = stats::rnorm(60))
  d$y <- 0.5 + 0.8 * d$x + stats::rnorm(60, 0, 0.7)
  fit <- frm(y ~ x, data = d,
             prior = set_prior("normal(0, 0.3)", class = "b"))
  b <- fit$estimates$beta
  sig <- exp(fit$estimates$betad[[1L]])
  data_ll <- sum(stats::dnorm(d$y, b[[1L]] + b[[2L]] * d$x, sig, log = TRUE))
  log_prior <- stats::dnorm(b[[2L]], 0, 0.3, log = TRUE)
  # the prior moved the slope, so the two values are far apart
  expect_gt(abs(log_prior), 1)
  expect_lt(abs(as.numeric(logLik(fit)) - data_ll) / abs(data_ll),
            sqrt(.Machine$double.eps))
})

test_that("the split equals the no-prior objective at the MAP estimates", {
  d <- lp_data()
  routes <- list(
    laplace = list(),
    profile = list(control = frmtmb_control(profile = TRUE)))
  f0 <- frm(y ~ x + (1 | g), data = d)
  for (r in names(routes)) {
    fp <- do.call(frm, c(list(y ~ x + (1 | g), data = d, prior = lp_prior()),
                         routes[[r]]))
    # the Laplace data objective at (beta, theta) of the MAP fit; a
    # profiled tape would maximize beta again, so f0 is the plain one
    pv <- unlist(fp$estimates[unique(names(f0$obj$par))])
    data_ll <- -f0$obj$fn(unname(pv))
    ll <- as.numeric(logLik(fp))
    expect_lt(abs(ll - data_ll) / abs(data_ll), sqrt(.Machine$double.eps),
              label = r)
    pen <- as.numeric(logLik_penalized(fp))
    lpr <- attr(logLik_penalized(fp), "log_prior")
    # the old value differs by the whole log prior
    expect_lt(abs((pen - data_ll) - lpr) / abs(lpr),
              sqrt(.Machine$double.eps), label = r)
    expect_gt(abs(lpr) / abs(data_ll), 1e-3)
  }
})

test_that("the split holds under quadrature", {
  d <- lp_data()
  set.seed(5)
  d$yb <- stats::rbinom(nrow(d), 1, stats::plogis(0.3 * d$x +
                                                    stats::rnorm(12)[d$g]))
  f0 <- frm(yb ~ x + (1 | g), data = d, family = bernoulli(),
            quadrature = TRUE)
  fp <- frm(yb ~ x + (1 | g), data = d, family = bernoulli(),
            quadrature = TRUE, prior = lp_prior())
  data_ll <- -f0$obj$fn(fp$opt$par)
  expect_lt(abs(as.numeric(logLik(fp)) - data_ll) / abs(data_ll),
            sqrt(.Machine$double.eps))
})

test_that("REML refuses logLik() for a prior on the integrated beta", {
  d <- lp_data()
  # the prior on sd reads an outer parameter: the split is exact
  ps <- set_prior("normal(0, 0.3)", class = "sd")
  f0 <- frm(y ~ x + (1 | g), data = d, REML = TRUE)
  fs <- frm(y ~ x + (1 | g), data = d, REML = TRUE, prior = ps)
  data_ll <- -f0$obj$fn(fs$opt$par)
  expect_lt(abs(as.numeric(logLik(fs)) - data_ll) / abs(data_ll),
            sqrt(.Machine$double.eps))
  # a prior on b is inside the integral REML computes
  fb <- frm(y ~ x + (1 | g), data = d, REML = TRUE,
            prior = set_prior("normal(0, 0.2)", class = "b"))
  expect_error(logLik(fb), class = "frmtmb_loglik_inner_prior")
  expect_error(logLik(fb), "REML integrates", fixed = TRUE)
  pen <- logLik_penalized(fb)
  expect_identical(as.numeric(pen), -fb$opt$objective)
  expect_true(is.na(attr(pen, "log_prior")))
  # summary() still prints, with no logLik to give
  s <- summary(fb)
  expect_true(is.na(s$loglik))
  out <- utils::capture.output(print(s))
  expect_true(any(grepl("inside the REML integral", out, fixed = TRUE)))
  expect_true(any(grepl("logLik_penalized:", out, fixed = TRUE)))
})

test_that("AIC(), BIC() and extractAIC() refuse a prior with a density", {
  d <- lp_data()
  f0 <- frm(y ~ x + (1 | g), data = d)
  fp <- frm(y ~ x + (1 | g), data = d,
            prior = set_prior("normal(0, 100)", class = "b"))
  expect_error(AIC(fp), class = "frmtmb_prior_ic")
  expect_error(BIC(fp), class = "frmtmb_prior_ic")
  expect_error(stats::extractAIC(fp), class = "frmtmb_prior_ic")
  expect_error(stats::step(fp, trace = 0), class = "frmtmb_prior_ic")
  # wherever the fit with the prior stands in the call
  expect_error(AIC(f0, fp), "model 2 has a prior", fixed = TRUE)
  expect_error(BIC(fp, f0), "model 1 has a prior", fixed = TRUE)
  # with no prior, and with a prior of bounds alone, the criteria answer
  fb <- frm(y ~ x + (1 | g), data = d,
            prior = set_prior("", class = "b", coef = "x", lb = -10))
  ll0 <- logLik(f0)
  aic0 <- -2 * as.numeric(ll0) + 2 * attr(ll0, "df")
  expect_lt(abs(AIC(f0) - aic0) / abs(aic0), sqrt(.Machine$double.eps))
  expect_lt(abs(AIC(fb) - AIC(f0)) / abs(AIC(f0)), 1e-8)
  expect_lt(abs(BIC(fb) - BIC(f0)) / abs(BIC(f0)), 1e-8)
  tab <- AIC(f0, fb)
  expect_identical(rownames(tab), c("f0", "fb"))
  expect_identical(as.numeric(logLik_penalized(fb)), as.numeric(logLik(fb)))
  expect_identical(attr(logLik_penalized(fb), "log_prior"), 0)
})

test_that("criteria that dispatch elsewhere read NA, never a number", {
  # AIC(lm, fp) dispatches on lm, AIC(logLik(fp)) on logLik: neither
  # reaches the refusal, so logLik() itself carries df = NA
  d <- lp_data()
  fl <- stats::lm(y ~ x, data = d)
  fp <- frm(y ~ x + (1 | g), data = d,
            prior = set_prior("normal(0, 0.1)", class = "b"))
  ll <- logLik(fp)
  expect_true(is.na(attr(ll, "df")))
  expect_match(attr(ll, "df_why"), "not a full degree of freedom",
               fixed = TRUE)
  expect_true(is.finite(as.numeric(ll)))
  expect_output(print(ll), "df=NA", fixed = TRUE)
  a <- AIC(fl, fp)
  expect_true(is.finite(a$AIC[1L]) && is.na(a$AIC[2L]))
  b <- BIC(fl, fp)
  expect_true(is.finite(b$BIC[1L]) && is.na(b$BIC[2L]))
  expect_true(is.na(AIC(ll)))
  # with no prior, df is a count and the same three spellings answer
  f0 <- frm(y ~ x + (1 | g), data = d)
  expect_false(is.na(attr(logLik(f0), "df")))
  expect_true(is.finite(AIC(logLik(f0))))
  expect_true(all(is.finite(AIC(fl, f0)$AIC)))
})

test_that("a uniform prior on a coefficient is bounds, not a density", {
  # its log density is a constant inside the bounds, so the estimates are
  # the maximum likelihood ones there and the criteria answer
  d <- lp_data()
  f0 <- frm(y ~ x + (1 | g), data = d)
  fu <- frm(y ~ x + (1 | g), data = d,
            prior = set_prior("uniform(-5, 5)", class = "b",
                                lb = -5, ub = 5))
  expect_lt(abs(fixef(fu, flatten = TRUE)[["x"]] -
                  fixef(f0, flatten = TRUE)[["x"]]) /
              abs(fixef(f0, flatten = TRUE)[["x"]]), 1e-6)
  expect_false(is.na(attr(logLik(fu), "df")))
  expect_lt(abs(AIC(fu) - AIC(f0)) / abs(AIC(f0)), 1e-8)
  # the maximized value differs from logLik() by the constant log(1 / 10)
  lpr <- attr(logLik_penalized(fu), "log_prior")
  expect_lt(abs(lpr - log(1 / 10)) / abs(log(1 / 10)),
            sqrt(.Machine$double.eps))
  expect_identical(summary(fu)$algorithm, "ML")
  # on a standard deviation the density is read through log(sd), whose
  # Jacobian is not constant: that is a density, and it is refused
  fs <- frm(y ~ x + (1 | g), data = d,
            prior = set_prior("uniform(0, 10)", class = "sd", ub = 10))
  expect_error(AIC(fs), class = "frmtmb_prior_ic")
  expect_true(is.na(attr(logLik(fs), "df")))
})

test_that("the split holds through every chart's Jacobian", {
  # one row per change of scale the prior is read through, beyond the sd
  # of the blocks above: a simplex, an AR coefficient, the gp length
  # scale and sd, and the Cholesky chart of an LKJ prior
  split_rel <- function(f0, fp) {
    expect_identical(names(f0$obj$par), names(fp$opt$par))
    data_ll <- -f0$obj$fn(fp$opt$par)
    expect_gt(abs(attr(logLik_penalized(fp), "log_prior")), 1e-3)
    abs(as.numeric(logLik(fp)) - data_ll) / abs(data_ll)
  }
  set.seed(4105)
  dmo <- data.frame(xo = sample(1:4, 120, TRUE))
  dmo$y <- c(0, 0.2, 1.4, 1.6)[dmo$xo] + stats::rnorm(120, 0, 0.5)
  po <- set_prior("dirichlet(c(2, 2, 2))", class = "simo", coef = "moxo1")
  expect_lt(split_rel(frm(y ~ mo(xo), data = dmo),
                      frm(y ~ mo(xo), data = dmo, prior = po)),
            sqrt(.Machine$double.eps), label = "mo()")
  set.seed(4103)
  dt <- data.frame(t = seq_len(80))
  dt$y <- as.numeric(stats::arima.sim(list(ar = 0.6), 80))
  pa <- set_prior("normal(0, 0.3)", class = "ar")
  expect_lt(split_rel(frm(y ~ 1 + ar(time = t), data = dt),
                      frm(y ~ 1 + ar(time = t), data = dt, prior = pa)),
            sqrt(.Machine$double.eps), label = "ar")
  set.seed(4102)
  dg <- data.frame(x = seq(0, 1, length.out = 30))
  dg$y <- sin(4 * dg$x) + stats::rnorm(30, 0, 0.3)
  pg <- set_prior("normal(0, 1)", class = "sdgp") +
    set_prior("inv_gamma(3, 0.5)", class = "lscale")
  expect_lt(split_rel(frm(y ~ gp(x), data = dg),
                      frm(y ~ gp(x), data = dg, prior = pg)),
            sqrt(.Machine$double.eps), label = "exact gp")
  set.seed(4107)
  dl <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, each = 10)))
  u <- matrix(stats::rnorm(32, 0, 0.6), 16)
  dl$y <- 1 + u[dl$g, 1] + (0.5 + u[dl$g, 2]) * dl$x + stats::rnorm(160)
  pl <- set_prior("lkj(4)", class = "cor") +
    set_prior("normal(0, 1)", class = "sd")
  expect_lt(split_rel(frm(y ~ x + (1 + x | g), data = dl),
                      frm(y ~ x + (1 + x | g), data = dl, prior = pl)),
            sqrt(.Machine$double.eps), label = "lkj")
})

test_that("logLik_penalized() on frm_multiple() refuses, as logLik() does", {
  set.seed(3)
  imps <- lapply(1:2, function(i) {
    d <- data.frame(x = stats::rnorm(40))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(40)
    d
  })
  m <- frm_multiple(bf(y ~ x), data = imps)
  expect_error(logLik_penalized(m), "no pooled answer", fixed = TRUE)
  expect_error(logLik(m), "no pooled answer", fixed = TRUE)
})

test_that("the nonlinear restart message names the penalized value", {
  # a and b start at the prior locations, 0, so the body a * b * x is 0
  # on every row: a saddle, and the restart reports what it compared
  set.seed(501)
  d <- data.frame(x = stats::runif(80, 0, 3))
  d$y <- 0.8 * d$x + stats::rnorm(80, 0, 0.3)
  msgs <- character(0)
  fit <- withCallingHandlers(
    frm(bf(y ~ a * b * x, a + b ~ 1, nl = TRUE), data = d,
        prior = set_prior("normal(0, 1)", nlpar = "a") +
          set_prior("normal(0, 1)", nlpar = "b")),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  hit <- grep("restarted from moved coefficients", msgs, value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, "penalized log-likelihood", fixed = TRUE)
  # and the number it prints last is the returned fit's maximized value,
  # not its logLik(), which is the log prior away (2.67 at this seed)
  got <- regmatches(hit, gregexpr("log-likelihood -?[0-9.]+", hit))[[1L]]
  nums <- as.numeric(sub("log-likelihood ", "", got, fixed = TRUE))
  pen <- as.numeric(logLik_penalized(fit))
  expect_lt(abs(nums[length(nums)] - pen) / abs(pen), 1e-5)
  expect_gt(abs(nums[length(nums)] - as.numeric(logLik(fit))) / abs(pen),
            1e-3)
})

test_that("summary() under a prior prints the log prior and no criteria", {
  d <- lp_data()
  fp <- frm(y ~ x + (1 | g), data = d, prior = lp_prior())
  s <- summary(fp)
  expect_identical(s$algorithm, "ML (MAP)")
  expect_true(is.na(s$AIC) && is.na(s$BIC))
  expect_identical(s$log_prior, attr(logLik_penalized(fp), "log_prior"))
  out <- utils::capture.output(print(s))
  expect_true(any(grepl("log prior:", out, fixed = TRUE)))
  expect_false(any(grepl("AIC:", out, fixed = TRUE)))
  # a prior of bounds alone is not MAP, and keeps the criteria
  fb <- frm(y ~ x + (1 | g), data = d,
            prior = set_prior("", class = "b", coef = "x", lb = -10))
  sb <- summary(fb)
  expect_identical(sb$algorithm, "ML")
  expect_false(is.na(sb$AIC))
})

test_that("the discrete check runs under a prior", {
  set.seed(1)
  d <- data.frame(y = stats::rpois(50, 3))
  # a mass that forgot a constant, so its log-likelihood is above 0
  bad <- frmtmb_family(
    "bad_pois", dpars = "mu", links = list(mu = "log"), type = "discrete",
    lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
      y * log(dpars$mu) - dpars$mu - lgamma(y + 1) + 10
    },
    init_dpars = list(mu = function(y, aterms) mean(y)))
  pr <- set_prior("normal(0, 0.5)", class = "Intercept")
  expect_warning(fit <- frm(y ~ 1, family = bad, data = d, prior = pr),
                 class = "frmtmb_discrete_loglik")
  expect_gt(as.numeric(logLik(fit)), 0)
  # the correct mass under the same prior does not warn
  expect_no_warning(frm(y ~ 1, family = poisson(), data = d, prior = pr))
  # and the negative control the old skip guarded: 50 zero counts and a
  # tight prior put the PENALIZED value above 0 (+3.35) while the data
  # value is below it (-0.34); a check on the penalized value warns here
  dz <- data.frame(y = rep(0L, 50))
  pz <- set_prior("normal(-5, 0.01)", class = "Intercept")
  expect_no_warning(fz <- frm(y ~ 1, family = poisson(), data = dz,
                              prior = pz))
  expect_gt(as.numeric(logLik_penalized(fz)), 0)
  expect_lt(as.numeric(logLik(fz)), 0)
})

test_that("the pooled likelihood rules refuse a prior with a density", {
  set.seed(3)
  imps <- lapply(1:2, function(i) {
    d <- data.frame(x = stats::rnorm(60), z = stats::rnorm(60))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(60)
    d
  })
  pr <- set_prior("normal(0, 1)", class = "b")
  m1 <- frm_multiple(bf(y ~ x + z), data = imps, prior = pr)
  m0 <- frm_multiple(bf(y ~ x), data = imps, prior = pr)
  expect_error(anova(m1, m0), "a pooled fit has a prior", fixed = TRUE)
  expect_error(anova(m1, m0, method = "D2"), "a pooled fit has a prior",
               fixed = TRUE)
  # a Wald rule reads no likelihood
  expect_s3_class(anova(m1, m0, method = "D1"), "data.frame")
  # and without the prior the likelihood rule answers
  n1 <- frm_multiple(bf(y ~ x + z), data = imps)
  n0 <- frm_multiple(bf(y ~ x), data = imps)
  expect_s3_class(anova(n1, n0), "data.frame")
})

test_that("frm_allfit() compares the value each optimizer maximized", {
  d <- lp_data()
  fp <- frm(y ~ x + (1 | g), data = d, prior = lp_prior())
  af <- frm_allfit(fp, optimizers = list(nlminb = "nlminb", optim = "optim"))
  pen <- as.numeric(logLik_penalized(fp))
  expect_identical(af$original$logLik, pen)
  expect_true(af$penalized)
  out <- utils::capture.output(print(af))
  expect_true(any(grepl("logLik_penalized()", out, fixed = TRUE)))
  # every optimizer lands where the original did, on the maximized value
  for (f in af$fits) {
    expect_lt(abs(as.numeric(logLik_penalized(f)) - pen) / abs(pen), 1e-6)
  }
  # with no prior the table is logLik() and says nothing more
  a0 <- frm_allfit(frm(y ~ x + (1 | g), data = d),
                   optimizers = list(nlminb = "nlminb"))
  expect_false(a0$penalized)
  expect_false(any(grepl("logLik_penalized()",
                         utils::capture.output(print(a0)), fixed = TRUE)))
})
