# Where a field's fit ends (R/gp-start.R): the mean's working marginal read
# again at the fit's optimum, the approximate field's grid there, and the
# warnings for a Laplace approximation that does not hold and for an
# approximate field whose basis is too coarse. Rows: dev/gpfield-rows.R,
# numbers in dev/gpfield-findings.md.

# lane arstep's rows (dev/gpstop-ng.R, seed 1)
gfd_ng <- function(case, n) {
  set.seed(74000 + 31 + n)
  x <- sort(round(stats::runif(n, 0, 10), 4))
  switch(case,
    pois_two = {
      eta <- sin(x) + 0.4 * sin(10 * x)
      data.frame(x = x, y = stats::rpois(n, exp(1.5 + eta / 2)))
    },
    student_step = data.frame(x = x, y = 2 * (x > 5) - 1 +
                                0.02 * stats::rt(n, 5)))
}

# lane gpstop's review sweep (dev/gpstop-fx2-sweep.R mk())
gfd_sweep <- function(body, fam, n, seed) {
  bodies <- c("exp_a", "exp_b", "logis_a", "mm_V", "mm_K", "pow_a", "pow_b")
  set.seed(81100 + 1000 * match(fam, c("gaussian", "poisson", "bernoulli",
                                       "student")) +
             100 * seed + n + 10 * match(body, bodies))
  x <- sort(round(stats::runif(n, 0, 10), 3))
  z <- round(stats::rnorm(n), 3)
  w <- round(stats::runif(n, 0.2, 5), 3)
  f <- if (seed %% 2) sin(x) else 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  eta <- switch(body,
    exp_a = (2 + fs) * exp(0.4 * z),
    logis_a = (3 + fs) / (1 + exp(-1.5 * z)),
    mm_V = (3 + fs) * w / (1 + w),
    pow_a = (2 + fs) * w^0.5)
  y <- switch(fam,
    gaussian = eta + 0.1 * stats::rnorm(n),
    poisson = stats::rpois(n, exp(0.6 * eta)),
    bernoulli = stats::rbinom(n, 1,
                              stats::plogis(1.5 * (eta - stats::median(eta)) /
                                              stats::sd(eta))))
  data.frame(x, z, w, y)
}

# the review's gate rows (dev/review-fitedge-designs.R rv_mk(, 2, 120, 2))
gfd_gate <- function(design) {
  di <- match(design, c("E", "F", "G", "H", "I", "J", "K", "L"))
  n <- 120L
  set.seed(52000 + 100 * 2 + n + 10 * di + 2)
  x <- sort(round(stats::runif(n, 0, 10), 3))
  z <- round(stats::rnorm(n), 3)
  g <- factor(sample(1:8, n, TRUE))
  fby <- factor(sample(c("p", "q"), n, TRUE))
  f <- 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  ge <- stats::rnorm(8, 0, 0.5)[g]
  y <- switch(design,
    E = (2 + fs) * exp(0.4 * z) + 0.1 * stats::rnorm(n),
    L = 2 * exp((0.3 + 0.2 * fs) * z) + 0.1 * stats::rnorm(n))
  data.frame(x, z, y)
}

test_that("a mean's field is read again at the fit's optimum", {
  # 0.82.0 stopped 1.59 below the regular maximum the start reaches, at code
  # 0 with no warning (backlog: the two rows of 0.80.0's Limits)
  skip_on_cran()
  d <- gfd_ng("pois_two", 80L)
  fit <- frm(bf(y ~ gp(x)), data = d, family = poisson())
  ref <- frm(bf(y ~ gp(x)), data = d, family = poisson(),
             start = list(theta = c(log(0.4), log(0.1))))
  l0 <- as.numeric(logLik(fit))
  l1 <- as.numeric(logLik(ref))
  expect_identical(ref$opt$convergence, 0L)
  expect_gte(l0, l1 - 1e-6 * abs(l1))
  expect_identical(fit$opt$convergence, 0L)
})

test_that("an approximate field gets a start at the fit's optimum", {
  # 0.82.0 stopped 1.18 below (b in the exponent) and at a white field 0.94
  # below a regular maximum (pow_a), each at code 0 with no warning
  skip_on_cran()
  d <- gfd_gate("L")
  fml <- bf(y ~ a * exp(b * z), a ~ 1, b ~ gp(x, k = 20), nl = TRUE)
  fit <- frm(fml, data = d, family = gaussian())
  ref <- frm(fml, data = d, family = gaussian(),
             start = list(theta = c(log(0.05), log(0.056))))
  l0 <- as.numeric(logLik(fit))
  l1 <- as.numeric(logLik(ref))
  expect_identical(ref$opt$convergence, 0L)
  expect_gte(l0, l1 - 1e-6 * abs(l1))
  d2 <- gfd_sweep("pow_a", "gaussian", 120L, 2L)
  fml2 <- bf(y ~ a * w^b, a ~ gp(x, k = 12), b ~ 1, nl = TRUE)
  fit2 <- frm(fml2, data = d2, family = gaussian())
  ref2 <- frm(fml2, data = d2, family = gaussian(),
              start = list(theta = c(log(0.8), log(0.9))))
  l2 <- as.numeric(logLik(ref2))
  expect_gte(as.numeric(logLik(fit2)), l2 - 1e-6 * abs(l2))
  # and it ends at a regular point, not at the white one
  bk <- fit2$frame$re_blocks[[1L]]
  expect_false(frmtmb:::gp_basis_white(bk, fit2$estimates$theta))
})

test_that("a field under its basis's spacing says so, and only there", {
  # the gaussian likelihood's maximum is there: k is too small for the
  # field, and only sdgp^2 * lscale is determined. Where the ridge is
  # proved level out to the white limit it is a boundary message, else
  # the warning (test-gp-white-limit.R reads which)
  d <- gfd_gate("E")
  said <- 0L
  fit <- allow_warnings(withCallingHandlers(
    frm(bf(y ~ a * exp(b * z), a ~ gp(x, k = 20), b ~ 1, nl = TRUE),
        data = d, family = gaussian()),
    frmtmb_gp_basis = function(cnd) {
      said <<- said + 1L
      if (inherits(cnd, "message")) invokeRestart("muffleMessage")
    }), "spacing its basis resolves")
  expect_identical(said, 1L)
  bk <- fit$frame$re_blocks[[1L]]
  expect_true(frmtmb:::gp_basis_white(bk, fit$estimates$theta))
  # the absent case: a regular approximate field, no warning
  d2 <- gfd_gate("L")
  fit2 <- frm(bf(y ~ a * exp(b * z), a ~ 1, b ~ gp(x, k = 20), nl = TRUE),
              data = d2, family = gaussian())
  expect_false(frmtmb:::gp_basis_white(fit2$frame$re_blocks[[1L]],
                                       fit2$estimates$theta))
})

test_that("the Laplace check warns where the approximation fails", {
  # student() on a step field: the first-order correction reads 2 to 3 nats
  # where importance sampling reads 2.72 (dev/gpfield-findings.md, 5125)
  skip_on_cran()
  # each row ends at a regular point or at a spike of the approximation,
  # by the rounding (st80 reaches the spike under OpenBLAS): either way it
  # warns, naming the inner Hessian at a spike and the correction elsewhere
  for (n in c(80L, 40L)) {
    d <- gfd_ng("student_step", n)
    msgs <- character()
    fit <- allow_warnings(
      withCallingHandlers(frm(bf(y ~ gp(x)), data = d, family = student()),
                          warning = function(w) {
                            msgs <<- c(msgs, conditionMessage(w))
                          }),
      c("The Laplace approximation", "Optimizer did not report", "gradient"),
      require = "The Laplace approximation")
    spike <- !frmtmb:::gp_inner_pd(fit$obj, fit$opt$par,
                                   frmtmb:::gp_relin_inner_ev_min)
    expect_true(any(grepl(if (spike) "inner Hessian is near singular" else
      "first-order correction", msgs, fixed = TRUE)), info = paste("n", n))
    if (!spike) {
      expect_gt(abs(frmtmb:::gp_laplace_tk(fit)), frmtmb:::gp_laplace_tk_max)
    }
  }
})

test_that("the Laplace check is quiet where the approximation holds", {
  # a white field on poisson counts is an observation-level effect, and
  # importance sampling agrees with the Laplace value there (-76.345 against
  # -76.324); a regular poisson field likewise
  d <- gfd_sweep("logis_a", "poisson", 40L, 1L)
  fit <- frm(bf(y ~ a / (1 + exp(-b * z)), a ~ gp(x), b ~ 1, nl = TRUE),
             data = d, family = poisson())
  tk <- frmtmb:::gp_laplace_tk(fit)
  expect_lt(abs(tk), frmtmb:::gp_laplace_tk_max)
  # the white field it is
  ls <- exp(fit$estimates$theta[2L])
  expect_lt(ls, 0.1 * stats::median(diff(sort(unique(d$x)))))
  skip_on_cran()
  d2 <- gfd_ng("pois_two", 80L)
  fit2 <- frm(bf(y ~ gp(x)), data = d2, family = poisson())
  expect_lt(abs(frmtmb:::gp_laplace_tk(fit2)), frmtmb:::gp_laplace_tk_max)
})

test_that("a white field on bernoulli rows warns", {
  # the Laplace value -63.97 stands about 40 above importance sampling's
  # there, where a regular maximum reads -67.07 by both (backlog: 0.81.0's
  # own first descent at a white-noise field)
  skip_on_cran()
  d <- gfd_sweep("exp_a", "bernoulli", 120L, 2L)
  allow_warnings(frm(bf(y ~ a * exp(b * z), a ~ gp(x), b ~ 1, nl = TRUE),
                     data = d, family = bernoulli()),
                 "Laplace approximation is not reliable",
                 require = "first-order correction")
})

test_that("a refit the sign probe returns is checked too", {
  # the probe's refit runs unannounced; the white field it returned on this
  # row read a Laplace value 22 above importance sampling's with no warning
  # (dev/reviews/2026-11-13-gpfield.md, F1)
  skip_on_cran()
  d <- gfd_sweep("mm_V", "bernoulli", 40L, 1L)
  fit <- allow_warnings(
    suppressMessages(frm(bf(y ~ V * w / (K + w), V ~ gp(x), K ~ 1, nl = TRUE),
                         data = d, family = bernoulli())),
    # under OpenBLAS 0.3.32 the refit also stops at code 8, a Newton step
    # from it gaining 1e-9, and says so
    c("Laplace approximation is not reliable", "Optimizer did not report"),
    require = "first-order correction")
  expect_gt(abs(frmtmb:::gp_laplace_tk(fit)), frmtmb:::gp_laplace_tk_max)
  # said once: a second ask of the same fit is quiet
  expect_no_warning(frmtmb:::gp_end_warn(fit))
})
