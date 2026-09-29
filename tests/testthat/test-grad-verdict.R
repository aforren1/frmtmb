# The convergence check on the gradient: what it must NOT warn on, and
# what it must still warn on.
#
# Every construction here was run against the 0.64.0 reference build
# first and the false-alarm blocks were SEEN TO FAIL there; the numbers
# are in dev/gradcheck-findings.md beside the script that produced them.
#
# No absolute tolerance appears below. Each assertion is a ratio to
# something the run measures: the fit's own `grad_tol`, its own standard
# errors, the bound it was given, this objective's own departure from a
# quadratic, `.Machine$double.eps`, or a quantity measured twice by two
# routes. Every bare number left is a MARGIN on one of those ratios.

grad_warned <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond))
    invokeRestart("muffleWarning")
  })
  list(fit = val, warnings = w,
       grad = any(grepl("Large maximum absolute gradient", w, fixed = TRUE)))
}

# --- a bound holding a parameter is not a convergence failure ---------

bounded_fit <- function(seed = 101, n = 300, cap = 0.1) {
  set.seed(seed)
  dd <- data.frame(x = rnorm(n))
  dd$y <- rnorm(n, 1 + 2 * dd$x, 1)
  # `cap` travels with the fit, so no block below repeats the number
  c(grad_warned(frm(bf(y ~ x), family = gaussian(), data = dd,
                    prior = set_prior("", class = "b", ub = cap))),
    list(cap = cap))
}

test_that("a parameter held by a bound does not raise the gradient warning", {
  r <- bounded_fit()
  # the behavioural assertion first, so a build without the new fields
  # fails on the WARNING rather than on a missing name
  expect_false(r$grad)
  d <- diagnose(r$fit, quiet = TRUE)
  tol <- r$fit$control$grad_tol
  # the constrained optimum IS the bound, so the old absolute criterion
  # was tripped by orders of magnitude
  expect_gt(d$max_grad / tol, 1e3)
  expect_identical(d$grad_bound_held, "x")
  # and the gradient over everything the bound does not hold is inside
  # the fit's own tolerance, so there is nothing to warn about
  expect_lt(d$grad_proj / tol, 1)
  # the parameter the verdict is ABOUT is not the one the raw gradient
  # names: that one is pinned by the bound and cannot move
  expect_identical(d$worst_grad, "x")
  expect_false(identical(d$grad_proj_par, d$worst_grad))
})

test_that("the bounded fit really is the constrained optimum", {
  r <- bounded_fit()
  nm <- frmtmb:::outer_par_names(r$fit)
  j <- match("x", nm)
  # ON the bound, to within the machine's own resolution at that value
  expect_lt(abs(unname(r$fit$opt$par[[j]]) - r$cap) / r$cap,
            sqrt(.Machine$double.eps))
  f0 <- as.numeric(r$fit$obj$fn(r$fit$opt$par))
  # every feasible move away from the bound is worse, in units of the
  # objective's own value at the optimum
  for (v in r$cap * c(0.99, 0.9, 0.5, 0)) {
    q <- r$fit$opt$par
    q[[j]] <- v
    expect_gt((as.numeric(r$fit$obj$fn(q)) - f0) / abs(f0), 0)
  }
})

test_that("diagnose() names the bound that holds the gradient", {
  r <- bounded_fit()
  out <- capture.output(diagnose(r$fit))
  expect_true(any(grepl("held by a bound", out, fixed = TRUE)))
  expect_true(any(grepl("the largest gradient no bound holds", out,
                        fixed = TRUE)))
  expect_true(any(grepl(diagnose(r$fit, quiet = TRUE)$grad_proj_par, out,
                        fixed = TRUE)))
})

test_that("the covariance report reads the full Hessian, not the free one", {
  # The defect the bounded fit runs into, stated as the implication that
  # WILL flip when the covariance becomes bound-aware, rather than as a
  # property of one cap. An earlier spelling asserted only
  # `expect_false(pdHess)`, which holds at cap 0.1 and 0.5 and fails at
  # every cap from 1 up, where the bound still holds `x` and the
  # machinery is still not bound-aware: it pinned how far the bound sits
  # from the unconstrained optimum, not bound-awareness
  # (dev/gradcheck-rev-10-flip.R). Filed in dev/test-backlog.md.
  r <- bounded_fit()
  d <- diagnose(r$fit, quiet = TRUE)
  expect_identical(d$grad_bound_held, "x")
  p <- r$fit$opt$par
  H <- r$fit$obj$he(p)
  H <- (H + t(H)) / 2
  free <- setdiff(seq_along(p),
                  match(d$grad_bound_held, frmtmb:::outer_par_names(r$fit)))
  ev_full <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  ev_free <- eigen(H[free, free, drop = FALSE], symmetric = TRUE,
                   only.values = TRUE)$values
  # the premise, measured on this fit rather than assumed: the full
  # Hessian is indefinite here and the free-set one is not. Skip rather
  # than fail if a future optimizer stops somewhere the premise fails,
  # because the implication below is only about the case it describes.
  skip_if_not(min(ev_full) < 0,
              "the full Hessian is positive definite at this optimum")
  expect_gt(min(ev_free) / max(ev_free), 0)
  # and the report follows the FULL one, which is the defect: flip this
  # to expect_true when the covariance is restricted to the free set
  expect_false(d$pdHess)
})

test_that("with the bound gone, nothing is excluded", {
  # The complement of the guard: no bound, so no component is held and
  # the projected gradient is the maximum absolute gradient IN THE UNITS
  # THE VERDICT USES. Those are the raw units only when the fit was not
  # standardized internally, because `max_grad` is raw while `grad_proj`
  # carries `par_units`; on an autoscaled fit they differ by the scaling
  # (1001 against 1.685 on one measured example, seed 9202). Asserting
  # `identical()` without checking that is asserting the fixture, so the
  # scaling is asserted first.
  set.seed(101)
  dd <- data.frame(x = rnorm(300))
  dd$y <- rnorm(300, 1 + 2 * dd$x, 1)
  f <- frm(bf(y ~ x), family = gaussian(), data = dd)
  d <- diagnose(f, quiet = TRUE)
  expect_identical(d$grad_bound_held, character(0))
  expect_null(f$par_units)
  expect_identical(d$grad_proj, d$max_grad)
  expect_identical(d$grad_proj_par, d$worst_grad)
})

test_that("grad_proj is in the units the verdict uses, max_grad is raw", {
  # the case the block above excludes, so the units claim is measured and
  # not assumed. A column spread far below one engages autoscale, which
  # gives the fit a par_units vector, and the two numbers then differ by
  # exactly that factor at the argmax.
  set.seed(9202)
  n <- 600
  dd <- data.frame(x = rnorm(n))
  dd$z <- dd$x * 1e-5
  dd$y <- rnorm(n, 1 + 2 * dd$x, 1)
  f <- frm(bf(y ~ z), family = gaussian(), data = dd,
           control = frmtmb_control(autoscale = TRUE))
  skip_if(is.null(f$par_units), "autoscale did not engage on this design")
  d <- diagnose(f, quiet = TRUE)
  expect_gt(max(f$par_units) / min(f$par_units), 1)
  # an IDENTITY, checked numerically: grad_proj is the raw gradient times
  # par_units, so with nothing held it is the max of that product
  gu <- abs(drop(f$obj$gr(f$opt$par)) * f$par_units)
  expect_identical(d$grad_proj, max(gu))
  expect_identical(d$max_grad, max(abs(drop(f$obj$gr(f$opt$par)))))
  # and the printed block warns the reader that the two lines differ
  expect_true(any(grepl("standardized internally",
                        capture.output(diagnose(f)), fixed = TRUE)))
})

test_that("an importance-corrected fit gets no verdict at all", {
  # check_convergence() refuses to judge a Monte Carlo gradient, so
  # diagnose() must not publish a headroom differenced from one
  set.seed(506)
  n <- 400
  dd <- data.frame(x = rnorm(n), g = factor(rep_len(1:40, n)))
  dd$y <- rpois(n, exp(0.3 + 0.6 * dd$x + rnorm(40, 0, 0.6)[dd$g]))
  f <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = poisson(),
                            data = dd, importance = 200L))
  skip_if(is.null(f$importance), "the importance correction did not run")
  d <- diagnose(f, quiet = TRUE)
  expect_true(is.finite(d$max_grad))
  expect_true(is.na(d$grad_proj))
  expect_true(is.na(d$grad_headroom))
  expect_identical(d$grad_bound_held, character(0))
  out <- capture.output(diagnose(f))
  expect_true(any(grepl("Monte Carlo gradient", out, fixed = TRUE)))
})

# The clean line on an importance fit, both ways round. The first spelling
# of the guard above withheld it unconditionally, which is a REGRESSION
# against 0.64.0: that build gated the line on the raw gradient, which a
# well-converged importance fit passes. The rule is that no verdict means
# the old gate decides, so this needs a fit on each side of `grad_tol` and
# asserts the two opposite outcomes. Seeds and designs from
# dev/gradcheck-rev-20-impclean.R.

imp_clean <- function(seed, ng, per, nd) {
  set.seed(seed)
  dd <- data.frame(g = factor(rep(seq_len(ng), per)))
  dd$x <- rnorm(nrow(dd))
  re <- rnorm(ng, 0, 0.5)
  dd$y <- rpois(nrow(dd), exp(0.3 + 0.5 * dd$x + re[dd$g]))
  f <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = poisson(),
                            data = dd, importance = nd, se = TRUE))
  list(fit = f, out = capture.output(suppressWarnings(diagnose(f))),
       d = suppressWarnings(diagnose(f, quiet = TRUE)))
}

test_that("an importance fit under grad_tol still prints the clean line", {
  r <- imp_clean(9503, ng = 20, per = 50, nd = 2000L)
  skip_if(is.null(r$fit$importance),
          "the importance correction did not run")
  tol <- r$fit$control$grad_tol
  # the premise, measured on this fit: nothing but the gradient could
  # withhold the line, and the gradient is inside the tolerance
  skip_if_not(r$d$max_grad < tol, "this draw did not land under grad_tol")
  expect_identical(r$fit$opt$convergence, 0L)
  expect_true(r$d$pdHess)
  expect_length(r$d$bad_se, 0L)
  expect_null(r$d$singular)
  expect_null(r$d$separation)
  expect_null(r$d$unbounded_dpar)
  expect_null(r$d$predictor_scale)
  expect_true(any(grepl("No convergence problems detected", r$out,
                        fixed = TRUE)))
})

test_that("an importance fit over grad_tol does not print the clean line", {
  # the ABSENT case of the same guard, so it cannot pass by printing the
  # line unconditionally either
  hit <- FALSE
  for (s in list(list(9501, 25, 40, 2000L), list(9502, 30, 30, 4000L))) {
    r <- imp_clean(s[[1L]], ng = s[[2L]], per = s[[3L]], nd = s[[4L]])
    if (is.null(r$fit$importance)) next
    if (r$d$max_grad <= r$fit$control$grad_tol) next
    hit <- TRUE
    expect_false(any(grepl("No convergence problems detected", r$out,
                           fixed = TRUE)))
  }
  skip_if_not(hit, "no draw landed above grad_tol")
})

test_that("a bound whose gradient points INWARD is not excluded", {
  # grad_bound_active() must read the sign, or a parameter that happens
  # to sit on a bound it is not held by would be dropped from the check
  par <- c(a = 0.1, b = 0.1)
  bd <- list(lower = c(-Inf, -Inf), upper = c(0.1, 0.1))
  # a positive gradient at an upper bound pushes the minimizer DOWN, off
  # the bound, so nothing holds it there
  expect_identical(unname(frmtmb:::grad_bound_active(par, c(-5, 5), bd)),
                   c(TRUE, FALSE))
  lo <- list(lower = c(0.1, 0.1), upper = c(Inf, Inf))
  expect_identical(unname(frmtmb:::grad_bound_active(par, c(5, -5), lo)),
                   c(TRUE, FALSE))
  # an open box holds nothing
  expect_identical(
    unname(frmtmb:::grad_bound_active(par, c(5, -5),
                                      list(lower = c(-Inf, -Inf),
                                           upper = c(Inf, Inf)))),
    c(FALSE, FALSE))
  # and a parameter well inside its bound is not "on" it
  expect_identical(
    unname(frmtmb:::grad_bound_active(c(a = 0, b = 0), c(-5, 5), bd)),
    c(FALSE, FALSE))
})

# --- a large sample does not make a correct fit warn ------------------

test_that("a gaussian GLM at n = 50000 trips the trip-wire and is clean", {
  set.seed(1040)
  n <- 50000
  dd <- data.frame(x = rnorm(n))
  dd$y <- rnorm(n, 1 + 2 * dd$x, 2)
  r <- grad_warned(frm(bf(y ~ x), family = gaussian(), data = dd))
  expect_false(r$grad)
  d <- diagnose(r$fit, quiet = TRUE)
  tol <- r$fit$control$grad_tol
  # the absolute gradient IS above the tolerance, so this fit reaches the
  # second stage rather than passing the first
  expect_gt(d$max_grad / tol, 1)
  # what is left on the table is a vanishing fraction of the tolerance
  expect_lt(d$grad_headroom / tol, 1)
  # and the fit is the closed form: every coefficient is inside a tiny
  # fraction of its own standard error of lm()'s
  lf <- lm(y ~ x, data = dd)
  se <- summary(lf)$coefficients[, 2]
  expect_lt(max(abs(coef(r$fit)[1:2] - coef(lf)) / se), 1e-3)
})

test_that("an ordinal fit at n = 6000 agrees with polr and stays clean", {
  skip_if_not_installed("MASS")
  set.seed(103)
  n <- 6000
  dd <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  eta <- 0.8 * dd$x1 - 0.5 * dd$x2
  dd$yo <- cut(eta + rlogis(n), breaks = c(-Inf, -1, 0.5, 2, Inf),
               labels = FALSE)
  dd$yof <- factor(dd$yo, ordered = TRUE)
  r <- grad_warned(frm(bf(yo ~ x1 + x2), family = cumulative(), data = dd))
  expect_false(r$grad)
  d <- diagnose(r$fit, quiet = TRUE)
  tol <- r$fit$control$grad_tol
  expect_gt(d$max_grad / tol, 1)
  expect_lt(d$grad_headroom / tol, 1)
  pl <- MASS::polr(yof ~ x1 + x2, data = dd, method = "logistic")
  # the two log likelihoods agree to far inside the tolerance the fit was
  # judged by, which is the point: the fit is at the same optimum an
  # independent implementation reaches
  gap <- abs(as.numeric(logLik(r$fit)) - as.numeric(logLik(pl)))
  expect_lt(gap / tol, 1)
})

# --- the headroom is the real remaining log likelihood ----------------

test_that("the headroom predicts what a Newton step actually gains", {
  set.seed(402)
  n <- 1500
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- rpois(n, exp(0.4 + 0.7 * dd$x - 0.4 * dd$z))
  r <- grad_warned(frm(
    bf(y ~ x + z), family = poisson(), data = dd,
    control = frmtmb_control(restarts = 0,
                             optCtrl = list(rel.tol = 1e-2, x.tol = 1e-2,
                                            iter.max = 1000,
                                            eval.max = 1000))))
  expect_true(r$grad)
  d <- diagnose(r$fit, quiet = TRUE)
  # measured twice by two routes: the criterion's prediction, and the
  # objective evaluated at the step it predicts
  p <- r$fit$opt$par
  g <- drop(r$fit$obj$gr(p))
  H <- r$fit$obj$he(p)
  ch <- chol((H + t(H)) / 2)
  step <- -backsolve(ch, backsolve(ch, g, transpose = TRUE))
  f0 <- as.numeric(r$fit$obj$fn(p))
  realized <- f0 - as.numeric(r$fit$obj$fn(p + step))
  half <- f0 - as.numeric(r$fit$obj$fn(p + 0.5 * step))
  # The yardstick is measured, not written down. Along the Newton
  # direction an exactly quadratic objective drops `D * (2t - t^2)`, so
  # the half-step drop is exactly 0.75 of the full one; the departure
  # from 0.75 is this objective's non-quadraticity AT THIS POINT, and the
  # headroom's own error is that same third-order term. Measured over six
  # constructions (dev/gradcheck-14-quadyard.R) the ratio of the two is
  # 0.96 to 1.29, so 10 is an order of magnitude of headroom. The floor
  # keeps an objective that is quadratic to the last bit from demanding
  # exactness of a floating-point difference.
  nonquad <- abs(half / realized / 0.75 - 1)
  yard <- max(nonquad, sqrt(.Machine$double.eps))
  expect_lt(abs(d$grad_headroom / realized - 1) / yard, 10)
})

# --- what must still warn --------------------------------------------

test_that("a loosened optimizer still warns, on a GLM and a GLMM", {
  set.seed(402)
  n <- 1500
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- rpois(n, exp(0.4 + 0.7 * dd$x - 0.4 * dd$z))
  loose <- frmtmb_control(restarts = 0,
                          optCtrl = list(rel.tol = 1e-2, x.tol = 1e-2,
                                         iter.max = 1000, eval.max = 1000))
  r1 <- grad_warned(frm(bf(y ~ x + z), family = poisson(), data = dd,
                        control = loose))
  # nlminb reports SUCCESS here, so the gradient check is the only thing
  # standing between the user and a fit short of its optimum
  expect_identical(r1$fit$opt$convergence, 0L)
  expect_true(r1$grad)

  set.seed(403)
  n <- 1200
  d3 <- data.frame(x = rnorm(n), g = factor(rep(1:40, 30)))
  d3$y <- rpois(n, exp(0.4 + 0.6 * d3$x + rnorm(40, 0, 0.7)[d3$g]))
  r2 <- grad_warned(frm(bf(y ~ x + (1 | g)), family = poisson(), data = d3,
                        control = loose))
  expect_identical(r2$fit$opt$convergence, 0L)
  expect_true(r2$grad)
  # both fits are short of their optimum by more than the tolerance, and
  # the check says so in log-likelihood units
  expect_gt(diagnose(r1$fit, quiet = TRUE)$grad_headroom /
              r1$fit$control$grad_tol, 1)
  expect_gt(diagnose(r2$fit, quiet = TRUE)$grad_headroom /
              r2$fit$control$grad_tol, 1)
})

test_that("a fit stopped at its iteration cap still warns", {
  set.seed(201)
  n <- 500
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- rnorm(n, 3 + 2 * dd$x - 1.5 * dd$z, 1)
  r <- grad_warned(frm(
    bf(y ~ x + z), family = gaussian(), data = dd,
    control = frmtmb_control(restarts = 0,
                             optCtrl = list(iter.max = 3, eval.max = 9))))
  expect_true(r$grad)
  ref <- frm(bf(y ~ x + z), family = gaussian(), data = dd)
  short <- as.numeric(logLik(ref)) - as.numeric(logLik(r$fit))
  expect_gt(short / r$fit$control$grad_tol, 1)
})

test_that("a start far away warns even when the curvature is unusable", {
  set.seed(201)
  n <- 500
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- rnorm(n, 3 + 2 * dd$x - 1.5 * dd$z, 1)
  r <- grad_warned(frm(
    bf(y ~ x + z), family = gaussian(), data = dd,
    start = list(beta = c(50, -40, 30)),
    control = frmtmb_control(restarts = 0,
                             optCtrl = list(iter.max = 1, eval.max = 3))))
  expect_true(r$grad)
  # the Hessian there is not positive definite, so the headroom cannot be
  # measured and the check warns rather than staying silent
  expect_true(any(grepl("curvature there is unusable", r$warnings,
                        fixed = TRUE)))
  expect_true(all(is.na(diagnose(r$fit, quiet = TRUE)$grad_headroom)))
})

test_that("a flat ridge warns on curvature the gradient does not show", {
  set.seed(502)
  n <- 1000
  x1 <- rnorm(n)
  dd <- data.frame(x1 = x1, x2 = x1 + rnorm(n, 0, 1e-5))
  dd$y <- rnorm(n, 1 + 2 * dd$x1, 1)
  r <- grad_warned(frm(bf(y ~ x1 + x2), family = gaussian(), data = dd,
                       control = frmtmb_control(restarts = 0)))
  expect_true(r$grad)
  # and it warns about the CURVATURE, not the gradient's size
  expect_true(any(grepl("in log-likelihood", r$warnings, fixed = TRUE)))
  d <- diagnose(r$fit, quiet = TRUE)
  tol <- r$fit$control$grad_tol
  # the gradient is barely over the trip-wire while the headroom is
  # hundreds of times the tolerance: the ill-conditioned direction is
  # what the second stage measures and the first cannot see
  expect_lt(d$max_grad / tol, 100)
  expect_gt(d$grad_headroom / tol, 100)
})

# --- grad_tol still drives both readings ------------------------------

test_that("grad_tol tightens and loosens the whole verdict", {
  set.seed(1030)
  n <- 1000
  dd <- data.frame(x = rnorm(n))
  dd$y <- rnorm(n, 1 + 0.8 * dd$x, 1)
  clean <- grad_warned(frm(bf(y ~ x), family = gaussian(), data = dd))
  expect_false(clean$grad)
  # a tolerance under the fit's own headroom brings the warning back
  tight <- grad_warned(frm(bf(y ~ x), family = gaussian(), data = dd,
                           control = frmtmb_control(grad_tol = 1e-14)))
  expect_true(tight$grad)

  set.seed(402)
  n <- 1500
  d2 <- data.frame(x = rnorm(n), z = rnorm(n))
  d2$y <- rpois(n, exp(0.4 + 0.7 * d2$x - 0.4 * d2$z))
  loose <- grad_warned(frm(
    bf(y ~ x + z), family = poisson(), data = d2,
    control = frmtmb_control(restarts = 0, grad_tol = 1e4,
                             optCtrl = list(rel.tol = 1e-2, x.tol = 1e-2,
                                            iter.max = 1000,
                                            eval.max = 1000))))
  # the same fit that warns at the default is silent at a tolerance above
  # its own gradient
  expect_false(loose$grad)
})
