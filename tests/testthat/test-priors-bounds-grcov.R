# Priors and bounds for sampling/fitting, and gr(cov=) known-covariance
# random effects.

sim_lmm <- function(seed = 301, n = 150, ng = 15) {
  set.seed(seed)
  dd <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng),
                                                length.out = n)))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(ng, 0, 0.7)[dd$g], 1)
  dd
}


test_that("prior name resolution: classes, coefficients, errors", {
  dd <- sim_lmm()
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  e <- frmtmb:::resolve_priors(fit, list(beta = prior_normal(0, 5)))
  expect_identical(e[[1]]$comp, "beta")
  expect_length(e[[1]]$idx, 2)
  expect_error(frmtmb:::resolve_priors(fit, list(zzz = prior_normal())),
               "Unknown parameter")
  expect_error(frmtmb:::resolve_priors(fit, list(x = 5)),
               "prior object")
})


test_that("hard bounds constrain the ML fit", {
  dd <- sim_lmm(seed = 303)
  fit0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  est0 <- fixef(fit0)$mu[["x"]]
  expect_lt(est0, 1)
  fitb <- suppressWarnings(
    frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
        prior = set_prior("", class = "b", coef = "x", lb = 1))
  )
  expect_equal(fixef(fitb)$mu[["x"]], 1, tolerance = 1e-6)
  expect_lt(as.numeric(logLik(fitb)), as.numeric(logLik(fit0)))
  expect_error(frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
                   prior = set_prior("", class = "b", coef = "zzz",
                                     lb = 0)),
               "Prior target not found")
})

test_that("gr(cov=) matches a hand-rolled correlated-intercepts reference", {
  set.seed(304)
  ng <- 30
  # a random PSD correlation matrix over the levels
  R <- crossprod(matrix(rnorm(ng * ng), ng)) / ng
  A <- stats::cov2cor(R)
  dimnames(A) <- list(as.character(seq_len(ng)),
                      as.character(seq_len(ng)))
  b_true <- drop(crossprod(chol(A), rnorm(ng))) * 0.8
  n <- 600
  g <- factor(rep(seq_len(ng), each = n / ng))
  x <- rnorm(n)
  dd <- data.frame(y = 1 + 0.5 * x + b_true[g] + rnorm(n, 0, 0.6),
                   x = x, g = g)

  fit <- frm(bf(y ~ x + (1 | gr(g, cov = A))) + gaussian(), data = dd)
  expect_identical(fit$frame$re_blocks[[1]]$covstruct, "gr_cov")

  yv <- dd$y; xv <- dd$x; gi <- as.integer(dd$g)
  Ad <- unname(A)
  nll_ref <- function(p) {
    Sigma <- exp(2 * p$lsd) * Ad
    nll <- -sum(RTMB::dmvnorm(p$u, 0, Sigma, log = TRUE))
    mu <- p$b[1] + p$b[2] * xv + p$u[gi]
    nll - sum(RTMB::dnorm(yv, mu, exp(p$ls), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(b = c(0, 0), ls = 0, lsd = 0,
                              u = numeric(ng)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)

  # sd recovered, methods work
  sd_hat <- sqrt(VarCorr(fit)[[1]][1, 1])
  expect_lt(abs(sd_hat - 0.8), 0.4)
  expect_identical(dim(ranef(fit)[[1]]), c(30L, 1L))
  cv <- confint_varcorr(fit)
  expect_true(all(cv$lwr < cv$estimate & cv$estimate < cv$upr))
  expect_equal(predict(fit, newdata = dd), predict(fit),
               tolerance = 1e-8)

  # validations
  A2 <- A; dimnames(A2) <- NULL
  expect_error(frm(bf(y ~ x + (1 | gr(g, cov = A2))) + gaussian(),
                   data = dd), "dimnames")
})

# ---------------------------------------------------------------------
# Prior-carried bounds address the same parameters the distribution path
# does. Before v0.49 a bound was keyed by the design-matrix COLUMN name,
# so every sub-formula spelled its intercept "(Intercept)": an
# nlpar/dpar/resp-addressed bound named no outer parameter and, with
# several such parameters, all of them collided on that one key.
# ---------------------------------------------------------------------

sim_nl_psy <- function(seed = 11, n = 400, guess = 0.25, thr = 5) {
  set.seed(seed)
  x <- runif(n, 0, 10)
  p <- guess + (1 - guess) * stats::plogis(x - thr)
  data.frame(y = stats::rbinom(n, 1, p), x = x, w = stats::rnorm(n))
}

psy_form <- function(rhs = "1") {
  bf(y ~ guess + (1 - guess) * plogis(x - thr), guess ~ 1,
     stats::as.formula(paste("thr ~", rhs)), nl = TRUE)
}

test_that("a prior bound on an nlpar lands on that parameter, and fits", {
  dd <- sim_nl_psy()
  form <- psy_form()
  st <- list(beta = c(0.5, 4))

  # an identity-link Bernoulli mean leaves [0, 1] on the way to the
  # optimum, so nlminb reports the evaluation it could not take
  fit <- suppressWarnings(
    frm(form, family = bernoulli(link = "identity"), data = dd,
        start = st,
        prior = set_prior("", nlpar = "guess", lb = 0, ub = 1)))
  # the defect refused this outright ("Unknown parameter(s) in bounds:
  # (Intercept)"), so merely fitting is half the regression
  expect_s3_class(fit, "frmtmb_fit")
  expect_lt(abs(unname(fixef(fit)$guess) - 0.25), 0.1)
  expect_lt(abs(unname(fixef(fit)$thr) - 5), 0.5)
  # and inside the box it was given, which is the point of the bound
  expect_gte(unname(fixef(fit)$guess), 0)
  expect_lte(unname(fixef(fit)$guess), 1)

  ri <- frmtmb:::resolve_prior_input(fit,
    set_prior("", nlpar = "guess", lb = 0, ub = 1))
  expect_identical(names(ri$lower), "guess_(Intercept)")
  expect_identical(unname(ri$lower), 0)
  expect_identical(unname(ri$upper), 1)

  # and it is the SAME box the argument spelling builds, which is what
  # lets the argument be the escape hatch rather than a second story
  nm <- outer_par_names(fit)
  expect_identical(frmtmb:::resolve_bounds(fit, ri$lower, ri$upper),
                   frmtmb:::resolve_bounds(fit, c(guess = 0), c(guess = 1)))
  expect_identical(nm[1], "guess_(Intercept)")
})

test_that("prior bounds on two nlpars do not collide on one key", {
  dd <- sim_nl_psy()
  fit <- suppressWarnings(
    frm(psy_form(), family = bernoulli(link = "identity"), data = dd,
        start = list(beta = c(0.5, 4))))
  ri <- frmtmb:::resolve_prior_input(fit,
    set_prior("", nlpar = "guess", lb = 0, ub = 1) +
      set_prior("", nlpar = "thr", lb = 2, ub = 8))
  expect_identical(ri$lower[["guess_(Intercept)"]], 0)
  expect_identical(ri$lower[["thr_(Intercept)"]], 2)
  expect_identical(ri$upper[["guess_(Intercept)"]], 1)
  expect_identical(ri$upper[["thr_(Intercept)"]], 8)
})

test_that("an nlpar bound covers every coefficient; coef= narrows it", {
  dd <- sim_nl_psy()
  fit <- suppressWarnings(
    frm(psy_form("1 + w"), family = bernoulli(link = "identity"),
        data = dd, start = list(beta = c(0.5, 4, 0))))
  # brms broadcasts a prior over class "b" of the named parameter, and a
  # bound is carried by the same specification, so it broadcasts too
  ri <- frmtmb:::resolve_prior_input(fit,
    set_prior("", nlpar = "thr", lb = -9, ub = 9))
  expect_setequal(names(ri$lower), c("thr_(Intercept)", "thr_w"))
  expect_true(all(ri$lower == -9))

  one <- frmtmb:::resolve_prior_input(fit,
    set_prior("", nlpar = "thr", coef = "w", lb = -9))
  expect_identical(names(one$lower), "thr_w")

  # the bare-name spelling of lower= stays a ONE-parameter alias and
  # still refuses the ambiguity that the prior vocabulary broadcasts over
  expect_error(frmtmb:::resolve_bounds(fit, c(thr = 0), NULL),
               "more than one coefficient")
})

test_that("prior bounds address dpar and resp coefficients", {
  set.seed(21)
  n <- 300
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x, exp(0.2 + 0.1 * dd$z))
  dd$y2 <- stats::rnorm(n, 0.3 - 0.4 * dd$x, 1)

  fd <- frm(bf(y ~ x, sigma ~ z) + gaussian(), data = dd)
  rid <- frmtmb:::resolve_prior_input(fd,
    set_prior("", dpar = "sigma", lb = -1, ub = 1))
  expect_identical(names(rid$lower), "sigma_z")
  # class "b" excludes the intercept, exactly as it does without a dpar
  rii <- frmtmb:::resolve_prior_input(fd,
    set_prior("", class = "Intercept", dpar = "sigma", lb = -0.5))
  expect_identical(names(rii$lower), "sigma_(Intercept)")
  expect_true(all(!is.na(match(names(rid$lower), outer_par_names(fd)))))

  fm <- frm(mvbf(bf(y ~ x) + gaussian(), bf(y2 ~ x) + gaussian()),
            data = dd)
  rim <- frmtmb:::resolve_prior_input(fm,
    set_prior("", resp = "y2", lb = -2, ub = 2))
  expect_identical(names(rim$lower), "y2_x")
  expect_true(all(!is.na(match(names(rim$lower), outer_par_names(fm)))))
})

test_that("class theta bounds one element or all of them", {
  dd <- sim_lmm(seed = 305)
  fit <- frm(bf(y ~ x + (1 + x | g)) + gaussian(), data = dd)
  nth <- length(fit$frame$par_template$theta)
  expect_gt(nth, 1L)

  all_th <- frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "theta", lb = -2, ub = 2))
  expect_length(all_th$lower, nth)

  # a single position, addressed by its outer name. This is the
  # per-element addressing that keeps class "theta" as expressive as the
  # argument spelling for the covariance parameters
  one <- frmtmb:::resolve_prior_input(fit,
    set_prior("", class = "theta", coef = "theta_2", lb = -1, ub = 1))
  expect_identical(names(one$lower), "theta_2")
  expect_identical(unname(one$lower), -1)
  expect_identical(unname(one$upper), 1)
})

test_that("a later bounds-only prior tightens an earlier one", {
  dd <- sim_nl_psy()
  form <- psy_form()
  st <- list(beta = c(0.5, 4))

  # "later wins" is the whole precedence rule now that frm() has no
  # lower/upper of its own. 0.4 is above the truth of 0.25, so the
  # tightened bound has to bite for the estimate to sit on it
  pr <- set_prior("", nlpar = "guess", lb = 0, ub = 1) +
    set_prior("", nlpar = "guess", lb = 0.4)
  fit <- suppressWarnings(
    frm(form, family = bernoulli(link = "identity"), data = dd,
        start = st, prior = pr))
  expect_equal(unname(fixef(fit)$guess), 0.4, tolerance = 1e-5)

  # the endpoint the later specification did not name is the earlier
  # one's, so tightening one side keeps the other
  ri <- frmtmb:::resolve_prior_input(fit, pr)
  expect_identical(unname(ri$lower["guess_(Intercept)"]), 0.4)
  expect_identical(unname(ri$upper["guess_(Intercept)"]), 1)
})

test_that("the retired lower=/upper= arguments are gone, not aliased", {
  dd <- sim_lmm(seed = 307)
  # a stale call must fail loudly rather than fit an unbounded model
  expect_error(
    frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, lower = c(x = 1)),
    "unused argument")
  expect_error(
    frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, upper = c(x = 1)),
    "unused argument")
  expect_false(any(c("lower", "upper") %in% names(formals(frm))))
})
