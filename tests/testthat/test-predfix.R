# Lane wt-predfix: five defects filed in the 0.62.0 round. Each block
# was seen to FAIL on the released 0.62.0 build (dev/predfix-log/).
# Records: dev/predfix-findings.md.

# --- fitted() on a multivariate fit ------------------------------------

pf_mv <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(20260921)
      n <- 150
      d <- data.frame(x = stats::rnorm(n))
      e <- matrix(stats::rnorm(2 * n), n, 2) %*%
        chol(matrix(c(1, 0.7, 0.7, 1), 2))
      d$y1 <- 1 + 0.5 * d$x + e[, 1]
      d$y2 <- -0.3 * d$x + e[, 2]
      cache <<- list(d = d, fit = frm(bf(mvbind(y1, y2) ~ x) + gaussian() +
                                        set_rescor(TRUE), data = d))
    }
    cache
  }
})

test_that("fitted() answers a multivariate fit in brms's n x 4 x nresp", {
  m <- pf_mv()
  out <- tryCatch(fitted(m$fit), error = function(e) conditionMessage(e))
  expect_true(is.array(out) && length(dim(out)) == 3L)
  # brms 2.23.0 on the same data: dim 150 x 4 x 2, rows unnamed, the
  # third dimension named by response (dev/predfix-log/brms.txt)
  expect_identical(dim(out), c(150L, 4L, 2L))
  expect_identical(dimnames(out),
                   list(NULL, c("Estimate", "Est.Error", "Q2.5", "Q97.5"),
                        c("y1", "y2")))
})

test_that("each layer of fitted(mv) is that response's own fitted()", {
  m <- pf_mv()
  out <- tryCatch(fitted(m$fit), error = function(e) NULL)
  for (r in c("y1", "y2")) {
    one <- fitted(m$fit, resp = r)
    expect_identical(if (is.null(out)) NULL else out[, , r], one, info = r)
  }
})

test_that("fitted(mv) takes brms's resp, dpar, scale and newdata", {
  m <- pf_mv()
  nd <- m$d[1:3, ]
  shape <- function(...) {
    tryCatch(dim(fitted(m$fit, ...)), error = function(e) conditionMessage(e))
  }
  # brms: one response is a matrix, several are an array in the order
  # asked for
  expect_identical(shape(resp = "y2"), c(150L, 4L))
  expect_identical(shape(resp = c("y1", "y2")), c(150L, 4L, 2L))
  expect_identical(shape(newdata = nd), c(3L, 4L, 2L))
  expect_identical(shape(dpar = "sigma"), c(150L, 4L, 2L))
  expect_identical(shape(scale = "linear"), c(150L, 4L, 2L))
  rev <- tryCatch(dimnames(fitted(m$fit, resp = c("y2", "y1")))[[3L]],
                  error = function(e) NULL)
  expect_identical(rev, c("y2", "y1"))
  err <- tryCatch(fitted(m$fit, resp = c("y1", "nope")),
                  error = function(e) e)
  expect_s3_class(err, "frmtmb_error")
})

test_that("fitted(mv) stacks a categorical response's layers as brms does", {
  set.seed(4)
  n <- 120
  d <- data.frame(x = stats::rnorm(n))
  d$y1 <- 1 + 0.5 * d$x + stats::rnorm(n)
  d$o <- cut(d$x + stats::rnorm(n), c(-Inf, -0.5, 0.5, Inf),
             labels = c("a", "b", "c"))
  f <- frm(mvbf(bf(y1 ~ x) + gaussian(), bf(o ~ x) + categorical()),
           data = d)
  out <- tryCatch(fitted(f), error = function(e) NULL)
  # brms 2.23.0 on these data (dev/predfix-log/brms-mvcat.txt): 120 x 4
  # x 4, layers y1, P(Y = a), P(Y = b), P(Y = c); resp reorders them
  expect_identical(dim(out), c(120L, 4L, 4L))
  expect_identical(dimnames(out)[[3L]],
                   c("y1", "P(Y = a)", "P(Y = b)", "P(Y = c)"))
  rev <- tryCatch(dimnames(fitted(f, resp = c("o", "y1")))[[3L]],
                  error = function(e) NULL)
  expect_identical(rev, c("P(Y = a)", "P(Y = b)", "P(Y = c)", "y1"))
  one <- fitted(f, resp = "o")
  expect_identical(dim(one), c(120L, 4L, 3L))
  if (!is.null(out)) {
    expect_identical(out[, , "P(Y = b)"], one[, , "P(Y = b)"])
    # brms's row 1 posterior means, 600 draws: P(Y = a) 0.1718,
    # P(Y = b) 0.5229, P(Y = c) 0.3053. ML against a posterior mean
    # agrees to the posterior's own spread (Est.Error about 0.05)
    brms_est <- c(0.1718, 0.5229, 0.3053)
    se <- out[1L, "Est.Error", 2:4]
    expect_true(all(abs(out[1L, "Estimate", 2:4] - brms_est) < se))
  }
})

# --- a quadrature fit's scalar Est.Error ----------------------------------

pf_quad <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(10)
      db <- data.frame(g = factor(rep(1:10, each = 6)), x = stats::rnorm(60))
      db$y <- stats::rbinom(60, 1, stats::plogis(db$x +
                                                   stats::rnorm(10)[db$g]))
      cache <<- list(d = db, fit = suppressWarnings(
        frm(bf(y ~ x + (1 | g)) + bernoulli(), data = db,
            quadrature = TRUE)))
    }
    cache
  }
})

test_that("fitted() on a quadrature fit says the group term is missing", {
  q <- pf_quad()
  nd <- q$d[1:2, c("x", "g")]
  # the scalar route, which lane wt-reunc found silent: 0.09407 against
  # 0.09340 at re_formula = NA, the group-effect term absent
  expect_warning(fitted(q$fit, newdata = nd),
                 class = "frmtmb_modes_conditional_se")
  expect_warning(fitted(q$fit),
                 class = "frmtmb_modes_conditional_se")
  expect_warning(fitted(q$fit, newdata = nd, scale = "linear"),
                 class = "frmtmb_modes_conditional_se")
})

test_that("fitted() on a quadrature fit is silent where nothing is missing", {
  q <- pf_quad()
  nd <- q$d[1:2, c("x", "g")]
  # a level the fit never saw carries the block variance, not a mode
  nd_new <- data.frame(x = nd$x, g = factor(c("new1", "new2")))
  expect_no_warning(fitted(q$fit, newdata = nd_new, allow_new_levels = TRUE))
  # re_formula = NA asks for no group effect, so none is left out
  expect_no_warning(fitted(q$fit, newdata = nd, re_formula = NA))
  # a Laplace fit carries the term
  fl <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = q$d)
  expect_no_warning(fitted(fl, newdata = nd))
})

# --- vcov() on a REML or profile fit with no random effects ---------------

pf_noreml <- function() {
  set.seed(5)
  d <- data.frame(x = stats::rnorm(80))
  d$cnt <- stats::rpois(80, exp(0.2 + 0.5 * d$x))
  d$y <- 1 + d$x + stats::rnorm(80)
  d
}

test_that("vcov() answers a poisson REML fit with no random effect", {
  d <- pf_noreml()
  f <- frm(bf(cnt ~ x) + poisson(), data = d, REML = TRUE)
  V <- tryCatch(vcov(f), error = function(e) conditionMessage(e))
  expect_true(is.matrix(V))
  # beta is the only parameter and it is integrated: its covariance is
  # the inverse of the inner Hessian, which is glm()'s once glm() is
  # converged as tightly (its default epsilon leaves 2e-5 in the ratio)
  Vg <- stats::vcov(stats::glm(cnt ~ x, family = stats::poisson(), data = d,
                              control = stats::glm.control(epsilon = 1e-12)))
  expect_equal(unname(if (is.matrix(V)) V else NA), unname(Vg),
               tolerance = 1e-5)
})

test_that("vcov() answers a poisson profile fit with no random effect", {
  d <- pf_noreml()
  f <- frm(bf(cnt ~ x) + poisson(), data = d,
           control = frmtmb_control(profile = TRUE))
  V <- tryCatch(vcov(f), error = function(e) conditionMessage(e))
  expect_true(is.matrix(V))
  Vg <- stats::vcov(stats::glm(cnt ~ x, family = stats::poisson(), data = d,
                              control = stats::glm.control(epsilon = 1e-12)))
  expect_equal(unname(if (is.matrix(V)) V else NA), unname(Vg),
               tolerance = 1e-5)
})

test_that("vcov() on a gaussian no-RE REML fit is lm()'s, as before", {
  d <- pf_noreml()
  f <- frm(bf(y ~ x) + gaussian(), data = d, REML = TRUE)
  expect_equal(unname(vcov(f)), unname(stats::vcov(stats::lm(y ~ x, d))),
               tolerance = 1e-5)
})

test_that("summary() and fitted() answer a poisson no-RE REML fit", {
  d <- pf_noreml()
  f <- frm(bf(cnt ~ x) + poisson(), data = d, REML = TRUE)
  s <- tryCatch(summary(f)$coefficients$mu[, "Std. Error"],
                error = function(e) conditionMessage(e))
  se_g <- sqrt(diag(stats::vcov(stats::glm(
    cnt ~ x, family = stats::poisson(), data = d,
    control = stats::glm.control(epsilon = 1e-12)))))
  expect_equal(unname(if (is.numeric(s)) s else NA), unname(se_g),
               tolerance = 1e-5)
  fe <- tryCatch(fitted(f)[, "Est.Error"], error = function(e) NULL)
  expect_true(is.numeric(fe) && all(is.finite(fe)))
})

# --- the default engages autoscale for a column spread below 1e-3 --------

pf_scaled <- function(s, seed = 23) {
  set.seed(seed)
  n <- 250
  xt <- stats::rnorm(n)
  data.frame(y = stats::rpois(n, exp(0.5 + 0.4 * xt)),
             g = 1 + 0.4 * xt + stats::rnorm(n), xs = xt * s)
}

test_that("a poisson y ~ 0 + x at scale 1e-6 reaches glm() by default", {
  d <- pf_scaled(1e-6)
  gl <- as.numeric(stats::logLik(
    stats::glm(y ~ 0 + xs, family = stats::poisson(), data = d)))
  f <- frm(bf(y ~ 0 + xs) + poisson(), data = d)
  # 0.62.0 stopped 62.578284 units short, convergence 0, no warning
  expect_equal(as.numeric(logLik(f)), gl, tolerance = 1e-8)
  expect_false(is.null(f$par_units))
})

test_that("a gaussian y ~ x at scale 1e-6 reaches lm() by default", {
  d <- pf_scaled(1e-6)
  ll <- as.numeric(stats::logLik(stats::lm(g ~ xs, data = d)))
  f <- frm(bf(g ~ xs) + gaussian(), data = d)
  expect_equal(as.numeric(logLik(f)), ll, tolerance = 1e-8)
})

test_that("the default leaves a fit with no column below 1e-3 bit for bit", {
  for (s in c(1, 1e-2, 1e6)) {
    d <- pf_scaled(s)
    f0 <- frm(bf(y ~ 0 + xs) + poisson(), data = d)
    f1 <- frm(bf(y ~ 0 + xs) + poisson(), data = d,
              control = frmtmb_control(autoscale = FALSE))
    expect_null(f0$par_units)
    expect_identical(f0$opt$par, f1$opt$par, info = format(s))
    expect_identical(logLik(f0), logLik(f1), info = format(s))
  }
})

test_that("the default leaves a REML fit's inner mu column alone", {
  # beta is integrated under REML = TRUE, and the inner Newton solver
  # was right at 1e-6 on every fit measured, so nothing engages
  d <- pf_scaled(1e-6)
  f0 <- frm(bf(g ~ xs) + gaussian(), data = d, REML = TRUE)
  f1 <- frm(bf(g ~ xs) + gaussian(), data = d, REML = TRUE,
            control = frmtmb_control(autoscale = FALSE))
  expect_null(f0$par_units)
  expect_identical(f0$opt$par, f1$opt$par)
})

test_that("autoscale = FALSE still turns the default off", {
  d <- pf_scaled(1e-6)
  f <- frm(bf(y ~ 0 + xs) + poisson(), data = d,
           control = frmtmb_control(autoscale = FALSE))
  expect_null(f$par_units)
  expect_null(frmtmb_control()$autoscale)
  expect_error(frmtmb_control(autoscale = "yes"), class = "frmtmb_error")
})

test_that("an engaged pre-fit warns no more than autoscale = FALSE", {
  # the pre-fit applies an Intercept bound to the CENTERED intercept,
  # where it binds, and warned on a fit that converged (punch round 1, B1)
  set.seed(1)
  n <- 200
  x <- stats::rnorm(n) * 1e-4 + 5e-4
  d <- data.frame(x = x, y = 1 + 2000 * (x - 5e-4) + stats::rnorm(n))
  count <- function(a, pr) {
    k <- 0L
    f <- withCallingHandlers(
      frm(bf(y ~ x), family = gaussian(), data = d, prior = pr,
          control = frmtmb_control(autoscale = a)),
      warning = function(w) {
        k <<- k + 1L
        invokeRestart("muffleWarning")
      })
    list(k = k, fit = f)
  }
  prs <- list(set_prior("", class = "Intercept", ub = 0.5),
              set_prior("", class = "b", coef = "x", lb = 3e4))
  for (i in seq_along(prs)) {
    pr <- prs[[i]]
    a <- count(NULL, pr)
    b <- count(FALSE, pr)
    # the Intercept bound engages; the lb bound's pre-fit stops at the
    # scaled bound without converging and the default falls back to the
    # FALSE fit (punch round 2), which is allowed and is what this pins
    if (i == 1L) expect_false(is.null(a$fit$par_units))
    expect_identical(a$k, b$k)
    # the reported fit is the one FALSE finds; a is at least as converged
    expect_equal(as.numeric(logLik(a$fit)), as.numeric(logLik(b$fit)),
                 tolerance = 1e-9)
  }
})

test_that("a random slope on a tiny column reaches the scale-1 optimum", {
  set.seed(2)
  g <- factor(rep(1:20, each = 15))
  x <- stats::rnorm(300)
  d <- data.frame(g = g, x0 = x,
                  y = 1 + 0.5 * x + stats::rnorm(20, 0, 0.7)[g] +
                    stats::rnorm(20, 0, 0.4)[g] * x + stats::rnorm(300))
  fo <- bf(y ~ x + (1 + x | g))
  d$x <- d$x0
  ref <- as.numeric(logLik(frm(fo, family = gaussian(), data = d)))
  # 0.62.0 and the first lane build: 1.88 units short at sd 1.07e-3
  # and 6.35 at 1.07e-6, code 0, because the slope's log sd started at
  # the scale-1 value and the X column alone was rescaled
  for (s in c(1e-2, 1e-3, 1e-6)) {
    d$x <- d$x0 * s
    f <- suppressWarnings(frm(fo, family = gaussian(), data = d))
    expect_equal(as.numeric(logLik(f)), ref, tolerance = 1e-8,
                 info = format(s))
    expect_false(is.null(f$par_units))
  }
})

test_that("verbose and diagnose() say what the default engaged", {
  d <- pf_scaled(1e-6)
  msgs <- character(0)
  f <- withCallingHandlers(
    frm(bf(y ~ 0 + xs) + poisson(), data = d,
        control = frmtmb_control(verbose = TRUE)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_true(any(grepl("^frmtmb: fit: poisson, ML, autoscale", msgs)))
  out <- utils::capture.output(diagnose(f))
  expect_false(any(grepl("refit with", out, fixed = TRUE)))
  expect_true(any(grepl("already standardized", out, fixed = TRUE)))
  f0 <- frm(bf(y ~ 0 + xs) + poisson(), data = d,
            control = frmtmb_control(autoscale = FALSE))
  out0 <- utils::capture.output(diagnose(f0))
  expect_true(any(grepl("autoscale = TRUE", out0, fixed = TRUE)))
})

# --- cs() in predict(): a silent wrong answer on 0.62.0 ------------------

pf_cs_data <- function() {
  set.seed(11)
  n <- 400
  x <- stats::rnorm(n)
  # the first threshold moves with x and the second against it, so the
  # effect is category-specific and not proportional
  p1 <- stats::plogis(-0.3 + 1.5 * x)
  p2 <- (1 - p1) * stats::plogis(0.8 - 1.2 * x)
  u <- stats::runif(n)
  data.frame(x = x, yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
}

# largest |z| of simulated category proportions against the exact
# probabilities, z in units of the binomial Monte Carlo error
pf_cs_maxz <- function(prop, p, nd) {
  max(abs(prop - p) / sqrt(pmax(p * (1 - p), 1e-12) / nd))
}

test_that("predict() carries cs() in sample and at newdata", {
  d <- pf_cs_data()
  nd <- data.frame(x = c(3, 3, -3, 0))
  ND <- 4000L
  for (fam in list(sratio(), acat(), cratio())) {
    f <- frm(bf(yo ~ cs(x)), family = fam, data = d)
    lab <- fam$family
    set.seed(1)
    pr <- predict(f, ndraws = ND, propagate_error = FALSE)
    # 0.62.0: max |z| 940.8 (sratio) and 533.4 (acat) in sample, every
    # row drawn as if cs(x) were absent (dev/predfix-log/cs-base.txt)
    expect_lt(pf_cs_maxz(pr, fitted(f)[, "Estimate", ], ND), 6,
              label = paste(lab, "in sample"))
    set.seed(1)
    pn <- predict(f, newdata = nd, ndraws = ND, propagate_error = FALSE)
    # 0.62.0: about (0.42, 0.39, 0.19) on every row, x = 3 and -3 alike
    expect_lt(pf_cs_maxz(pn, fitted(f, newdata = nd)[, "Estimate", ], ND),
              6, label = paste(lab, "at newdata"))
  }
})

# --- the default must never be less diagnostic than autoscale = FALSE ----

pf_count <- function(expr) {
  k <- 0L
  out <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(w) {
      k <<- k + 1L
      invokeRestart("muffleWarning")
    })
  list(value = out, warnings = k,
       error = inherits(out, "error"))
}

test_that("a separated fit warns under the default as under FALSE", {
  # the reviewer's design (punch round 2): the pre-fit ran to the
  # separated point, the reported fit started there and stopped with
  # code 0, no warning and Est.Error 2.8e132
  for (seed in 511:514) {
    set.seed(seed)
    d <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
    d$yb <- as.integer(d$z > 0)
    a <- pf_count(frm(yb ~ z + x, family = bernoulli(), data = d))
    b <- pf_count(frm(yb ~ z + x, family = bernoulli(), data = d,
                      control = frmtmb_control(autoscale = FALSE)))
    expect_false(a$error, info = seed)
    expect_gte(a$warnings, b$warnings)
    if (!a$error && !b$error) {
      # fallen back: the reported fit is the FALSE fit
      expect_identical(a$value$opt$par, b$value$opt$par)
    }
  }
})

test_that("a pre-fit that errors falls back under the default", {
  # seed 533 of the reviewer's random-slope design: the pre-fit died on
  # "NA/NaN gradient evaluation" and the default raised it, where
  # autoscale = FALSE returns a fit
  set.seed(533)
  g <- factor(rep(1:20, each = 12))
  d <- data.frame(g, x = stats::rnorm(240) * 0.03, z = stats::rnorm(240))
  d$yb <- as.integer(d$z > 0)
  fo <- yb ~ z + x + (1 + x | g)
  a <- pf_count(frm(fo, family = bernoulli(), data = d))
  b <- pf_count(frm(fo, family = bernoulli(), data = d,
                    control = frmtmb_control(autoscale = FALSE)))
  expect_false(a$error)
  expect_false(b$error)
  expect_gte(a$warnings, b$warnings)
  msgs <- character(0)
  withCallingHandlers(
    tryCatch(suppressWarnings(frm(fo, family = bernoulli(), data = d,
                                  control = frmtmb_control(verbose = TRUE))),
             error = function(e) NULL),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_true(any(grepl("autoscale pre-fit (failed|did not converge)",
                        msgs)))
  # asked for, the error stands and names the way out
  e <- tryCatch(frm(fo, family = bernoulli(), data = d,
                    control = frmtmb_control(autoscale = TRUE)),
                error = function(e) e, warning = function(w) w)
  if (inherits(e, "error")) {
    expect_match(conditionMessage(e), "autoscale = FALSE", fixed = TRUE)
  } else {
    expect_match(conditionMessage(e), "autoscale", fixed = TRUE)
  }
})

test_that("an engaged nonlinear fit announces its prior-placed start once", {
  # the pre-fit's template replaced the cold start and the message with
  # it (punch round 2, minor 3); autoscale = FALSE prints it once
  set.seed(7)
  n <- 200
  x <- stats::rnorm(n) * 1e-4
  t <- stats::runif(n, 0, 5)
  d <- data.frame(t, x, y = 3 * exp(-0.7 * t) + 2000 * x +
                    stats::rnorm(n, 0, 0.1))
  nlf <- bf(y ~ a * exp(-k * t) + c, a ~ 1, k ~ 1, c ~ 0 + x, nl = TRUE)
  pr <- set_prior("normal(0.5, 1)", nlpar = "k")
  count <- function(a) {
    m <- character(0)
    f <- withCallingHandlers(
      suppressWarnings(frm(nlf, data = d, prior = pr,
                           control = frmtmb_control(autoscale = a))),
      message = function(x) {
        m <<- c(m, conditionMessage(x))
        invokeRestart("muffleMessage")
      })
    list(n = sum(grepl("placed at the prior locations", m)),
         engaged = !is.null(f$par_units))
  }
  a <- count(NULL)
  b <- count(FALSE)
  expect_true(a$engaged)
  expect_identical(b$n, 1L)
  expect_identical(a$n, 1L)
})

test_that("a sound seeded fit is kept where the plain fit stalls", {
  # skew_normal alpha ~ 0 + xs at 1e-6: the pre-fit stops with a flat
  # direction, the plain fit stalls 12.09 units short with a gradient
  # warning, and the fit seeded from the pre-fit is a verified optimum
  # (punch round 2: refusing it lost 9.6 to 14.5 units on 4 seeds)
  mk <- function(s) {
    set.seed(4)
    xs <- -abs(stats::rnorm(250)) * 3
    yy <- xs + (abs(stats::rnorm(250)) - sqrt(2 / pi)) * 1.5 +
      stats::rnorm(250, 0, 0.3)
    data.frame(y = yy, xs = xs * s, xr = xs)
  }
  fo <- bf(y ~ xr, sigma ~ 1, alpha ~ 0 + xs)
  ref <- as.numeric(logLik(frm(fo, family = skew_normal(), data = mk(1))))
  a <- pf_count(frm(fo, family = skew_normal(), data = mk(1e-6)))
  expect_false(a$error)
  expect_equal(as.numeric(logLik(a$value)), ref, tolerance = 1e-8)
  expect_identical(a$warnings, 0L)
})
