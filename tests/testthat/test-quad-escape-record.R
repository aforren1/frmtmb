# A quadrature fit escapes a family's stationary point on the Laplace
# objective (quad_fit()) and again on the rule (fit_assembled()); the
# second escape's record replaced the first's, so a fit the Laplace escape
# moved and the rule's left alone reported a gain of 0
# (stationary_escape_keep(), release review of 0.82.0, minor 6).

test_that("the escape that moved a quadrature fit keeps its record", {
  lap <- c(starts = 2, gain = 24.27)
  # the rule's escape ran and gained nothing: the Laplace record stays
  opt <- list(objective = 1, stationary_escape = c(starts = 2, gain = 0))
  got <- stationary_escape_keep(opt, lap)
  expect_identical(got$stationary_escape, lap)
  expect_identical(got$stationary_escape_laplace, lap)
  # the rule's escape moved the fit: its record stays, the Laplace one
  # beside it
  rule <- c(starts = 2, gain = 0.5)
  opt$stationary_escape <- rule
  got <- stationary_escape_keep(opt, lap)
  expect_identical(got$stationary_escape, rule)
  expect_identical(got$stationary_escape_laplace, lap)
})

test_that("a restart that raises keeps its last finite-gradient point", {
  obj <- list(par = c(a = 0), env = new.env(),
              fn = function(x) -x[[1]],
              gr = function(x) if (x[[1]] > 5) NaN else -1)
  tr <- escape_track(obj)
  for (a in c(1, 3, 4, 7)) {
    tr$obj$fn(c(a = a))
    tr$obj$gr(c(a = a))
  }
  # a gradient at a point other than the last one read is not remembered
  tr$obj$fn(c(a = 2))
  tr$obj$gr(c(a = 4.5))
  b <- tr$best(simpleError("NA/NaN gradient evaluation"))
  expect_identical(b$par, c(a = 4))
  expect_identical(b$objective, -4)
  expect_identical(b$convergence, 1L)
  expect_match(b$message, "NA/NaN gradient evaluation", fixed = TRUE)
  # nothing remembered: the restart is dropped, as before
  expect_null(escape_track(obj)$best(simpleError("x")))
})

test_that("a restart back on the point does not veto an abandoned one", {
  # The restart from +2 reaches the maximum and then raises, as a run-off
  # does where the gradient stops being finite; the one from -2 comes back
  # to alpha = 0 with an objective 1e-12 lower, the sign of rounding
  # noise. A build that counted any lower objective as a gain dropped the
  # abandoned point and kept the stall (review of lane quadesc, finding 3,
  # dev/reviews/quadesc-review-gainhole.R). The abandoned result also
  # carries its evaluation and trial counts (finding 5)
  set.seed(4041)
  dl <- 4 / sqrt(17)
  dd <- data.frame(y = 1 + dl * abs(stats::rnorm(200)) +
                     sqrt(1 - dl^2) * stats::rnorm(200))
  real <- optimize_obj
  st <- new.env()
  alpha_of <- function(p) unname(p[names(p) == "betad"][2])
  local_mocked_bindings(optimize_obj = function(obj, control,
                                                bounds = list(lower = -Inf,
                                                              upper = Inf),
                                                par_units = NULL,
                                                verbose = 0L,
                                                start_par = obj$par,
                                                tally = NULL) {
    a <- alpha_of(start_par)
    r <- real(obj, control, bounds, par_units, verbose, start_par, tally)
    if (isTRUE(all.equal(a, 2))) stop("NA/NaN gradient evaluation")
    if (isTRUE(all.equal(a, -2))) {
      r <- st$stall
      r$objective <- r$objective - 1e-12
    }
    if (is.null(st$stall)) st$stall <- r
    r
  })
  fit <- allow_warnings(
    frm(bf(y ~ 1), family = skew_normal(), data = dd,
        start = list(betad = c(0, 0))),
    allowed = "Optimizer did not report convergence",
    require = "Optimizer did not report convergence")
  # the stall was on the point, and the escape left it
  expect_lt(abs(alpha_of(st$stall$par)), 0.05)
  expect_gt(abs(alpha_of(fit$opt$par)), 1)
  expect_identical(fit$opt$convergence, 1L)
  expect_match(fit$opt$message, "restart off a stationary point",
               fixed = TRUE)
  expect_gt(fit$opt$stationary_escape[["gain"]], 0)
  expect_false(is.null(fit$opt$evals))
  expect_false(is.null(fit$opt$nonfinite_trials))
})

test_that("the rule's escape takes alpha off its point on every rounding", {
  # Each restart of the rule's escape runs alpha far toward its end (the
  # half-normal limit), where the rule's gradient stops being finite.
  # Whether the optimizer stops first or raises is rounding's choice, and
  # a restart that raised was dropped: 0.83.0 stopped this fit (alpha ~
  # x, both coefficients started at 0) on alpha = 0 with the Laplace
  # escape's record on 4 of 21 roundings under the reference BLAS, 7 of
  # 21 under OpenBLAS 0.3.26 and 0.3.32, and on the ubuntu-26.04 runner
  # unscaled (item 5015, dev/quadesc-findings.md). The fit runs on the
  # data the runner saw (k = 0) and with the response scaled by 1 + 2e-15
  # * 9, a rounding on which 0.83.0 stopped there under all three local
  # builds. The data of the two designs before it in
  # dev/reviews/quadseam2-review-esc2.R (ALPHAX seed 11) are drawn first,
  # as there
  skip_on_cran()
  sn <- function(n, alpha) {
    dl <- alpha / sqrt(1 + alpha^2)
    dl * (abs(stats::rnorm(n)) - sqrt(2 / pi)) +
      sqrt(1 - dl^2) * stats::rnorm(n)
  }
  set.seed(1011)
  stats::rnorm(60 + 15 + 60)
  stats::rnorm(50 + 10)
  sn(50, 1.5)
  d0 <- data.frame(x = stats::rnorm(80), g = factor(rep(1:16, each = 5)))
  d0$y <- 0.5 * d0$x + stats::rnorm(16, 0, 1)[d0$g] + sn(80, 2)
  fo <- bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ x)
  # the exact marginal log-likelihood (one integral per group, by the
  # trapezoid rule over 20001 points spanning 12 group sds either way),
  # which reads neither the rule nor the Laplace approximation
  exact <- function(par, d) {
    b <- par[names(par) == "beta"]
    bd <- par[names(par) == "betad"]
    sdg <- exp(par[names(par) == "theta"][1])
    mu <- b[1] + b[2] * d$x
    al <- bd[2] + bd[3] * d$x
    u <- seq(-12 * sdg, 12 * sdg, length.out = 20001L)
    tot <- 0
    for (g in levels(d$g)) {
      l <- stats::dnorm(u, 0, sdg, log = TRUE)
      for (i in which(d$g == g)) {
        l <- l + RTMBdist::dskewnorm2(d$y[i], mu[i] + u, exp(bd[1]), al[i],
                                      log = TRUE)
      }
      w <- exp(l - max(l))
      tot <- tot + max(l) +
        log((u[2] - u[1]) * (sum(w) - 0.5 * (w[1] + w[length(w)])))
    }
    tot
  }
  for (k in c(0, 9)) {
    d <- d0
    d$y <- d$y * (1 + 2e-15 * k)
    # alpha stops far out on a likelihood that is nearly level there, and
    # each ending says it is not a maximum: the optimizer's code, a
    # gradient above grad_tol, a rule that re-anchoring cannot settle, or
    # standard errors that are not available (63 roundings over three
    # BLAS builds, dev/quadesc-findings.md, section 5). Any other warning
    # reaches the test
    fit <- allow_warnings(
      frm(fo, family = skew_normal(), data = d, quadrature = TRUE,
          start = list(betad = c(0, 0, 0))),
      allowed = c("Optimizer did not report convergence",
                  "Large maximum absolute gradient at the optimum",
                  "quadrature = TRUE did not settle the Gauss-Kronrod rule",
                  "Standard errors are not available for"))
    esc <- fit$opt[["stationary_escape"]]
    lap <- fit$opt[["stationary_escape_laplace"]]
    # the Laplace escape ran and its record is kept; the record the fit
    # reports is the rule's when the rule's escape gained, or that one
    expect_false(is.null(lap))
    expect_gt(lap[["gain"]], 0)
    expect_true(identical(esc, lap) || isTRUE(esc[["gain"]] > 0))
    # alpha runs toward its end
    a <- fit$estimates$betad[grep("^alpha", names(fit$estimates$betad))]
    expect_true(all(abs(a) > 1), info = paste("k =", k))
    # and the fit is not short against the Laplace fit, which escapes on
    # its own objective, both scored exactly
    lap_fit <- frm(fo, family = skew_normal(), data = d,
                   start = list(betad = c(0, 0, 0)))
    expect_gt(exact(fit$opt$par, d), exact(lap_fit$opt$par, d))
  }
})
