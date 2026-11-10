# The standard-error check on fits that carry a convergence warning
# (lane illse, item c). It used to skip them, so on a ridge whether a
# fit's standard errors were judged depended on which side of code 0
# rounding put nlminb (lane ci691). Now every fit is checked, and the
# warning says whether the stop is the maximum: a Newton step from it
# gains at most grad_tol, or more.

conv_conds <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

conv_phrase <- "Optimizer did not report convergence"

# nl `a + b`: the two enter only through their sum, an exact alias the
# check names on a converged fit
conv_alias_data <- function(seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(100))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(100, 0, 0.5)
  d
}
conv_alias_bf <- bf(y ~ a + b + c * x, a + b + c ~ 1, nl = TRUE)

test_that("a stop at the maximum nlminb did not call converged is judged", {
  d <- conv_alias_data()
  ref <- conv_conds(frm(conv_alias_bf, data = d))
  expect_equal(ref$value$opt$convergence, 0L)
  lost0 <- sdr_of(ref$value)$se_lost
  expect_true(length(lost0) > 0L)
  # the same model started at its own maximum, with no iteration
  # allowed: nlminb reports the iteration limit at the maximum
  st <- ref$value$estimates[c("beta", "betad")]
  r <- conv_conds(frm(conv_alias_bf, data = d, start = st,
                      control = frmtmb_control(optCtrl = list(
                        iter.max = 0))))
  expect_false(r$value$opt$convergence == 0L)
  w <- grep(conv_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "so the estimates are at the maximum", fixed = TRUE)
  # judged as the converged fit is: the same lost set, and the same
  # message for it (the nonlinear alias's own), beside the convergence
  # warning and no second one in other words
  expect_identical(sdr_of(r$value)$se_lost, lost0)
  expect_identical(setdiff(r$warnings, w), ref$warnings)
  expect_length(grep("Large maximum absolute gradient", r$warnings,
                     fixed = TRUE), 0L)
})

test_that("a stop short of the maximum says so and keeps its SEs", {
  # A curvature verdict at a point short of the maximum is not one about
  # the data, so such a fit keeps sdreport()'s standard errors, and the
  # warning says they are those of the stop
  # A poisson fit from an explicit start far from the maximum, so that
  # one iteration stops short whatever start frm() would choose itself
  # (from a least-squares start a gaussian y ~ x is at the maximum
  # after one iteration)
  set.seed(2)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.5 * d$x))
  r <- conv_conds(frm(y ~ x, data = d, family = poisson(),
                      start = list(beta = c(-3, 3)),
                      control = frmtmb_control(optCtrl = list(
                        iter.max = 1), restarts = 0)))
  expect_false(r$value$opt$convergence == 0L)
  w <- grep(conv_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(w, 1L)
  expect_false(grepl("so the estimates are at the maximum", w,
                     fixed = TRUE))
  expect_true(grepl("short of the maximum", w, fixed = TRUE) ||
                grepl("curvature at the stop is unusable", w, fixed = TRUE))
  # sdreport()'s covariance, not repaired or judged
  s <- sdr_of(r$value)
  expect_null(s$cov_fixed_prop)
  expect_length(s$se_lost, 0L)
  expect_true(all(is.finite(sqrt(diag(s$cov.fixed)))))
})

# An optimizer that runs nlminb and then reports its stop as nlminb's
# "false convergence (8)", so that a code-8 stop is made on every
# platform: at the maximum with `iter` large, short of it with `iter`
# small. A real code-8 stop is a rounding's choice (dev/code8-battery.R
# makes 109 of 319 fits stop there at the maximum)
conv_fake8 <- function(iter = 1000L) {
  function(par, fn, gr, lower, upper, control) {
    o <- stats::nlminb(par, fn, gr, lower = lower, upper = upper,
                       control = list(iter.max = iter, eval.max = 1000))
    list(par = o$par, objective = o$objective, convergence = 1L,
         message = "false convergence (8)")
  }
}

conv_classes <- function(expr) {
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

test_that("false convergence (8) at the maximum is a message of its class", {
  # Up to 0.80.0 it warned "Optimizer did not report convergence" on
  # every such fit, including the 109 code-8 stops of the bigx battery
  # that are all at the maximum (Newton gains 2.4e-11 to 9.9e-6)
  set.seed(2)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.5 * d$x))
  r <- conv_classes(frm(y ~ x, data = d, family = poisson(),
                        control = frmtmb_control(optimizer = conv_fake8())))
  expect_identical(r$value$opt$convergence, 1L)
  expect_length(r$warnings, 0L)
  at <- Filter(function(x) inherits(x, "frmtmb_converged_at_max"),
               r$messages)
  expect_length(at, 1L)
  expect_match(conditionMessage(at[[1L]]), "false convergence (8)",
               fixed = TRUE)
  expect_true(all(is.finite(sqrt(diag(vcov(r$value))))))
  # the guard's absent cases: a stop short of the maximum keeps the
  # warning (its class is frmtmb_conv_code) and gives no message
  s <- conv_classes(frm(y ~ x, data = d, family = poisson(),
                        start = list(beta = c(-3, 3)),
                        control = frmtmb_control(optimizer = conv_fake8(1L),
                                                 restarts = 0)))
  expect_true(any(vapply(s$warnings, inherits, NA, "frmtmb_conv_code")))
  expect_false(any(vapply(s$messages, inherits, NA,
                          "frmtmb_converged_at_max")))
  # at the maximum with a standard error lost (the alias a + b) it is
  # not known to be an interior maximum, and the warning stays
  da <- conv_alias_data()
  a <- conv_classes(frm(conv_alias_bf, data = da,
                        control = frmtmb_control(optimizer = conv_fake8())))
  expect_true(length(sdr_of(a$value)$se_lost) > 0L)
  expect_true(any(vapply(a$warnings, inherits, NA, "frmtmb_conv_code")))
  expect_false(any(vapply(a$messages, inherits, NA,
                          "frmtmb_converged_at_max")))
  # and so it does where the standard errors are not asked about
  q <- conv_classes(frm(y ~ x, data = d, family = poisson(),
                        control = frmtmb_control(optimizer = conv_fake8(),
                                                 check_se = "ignore")))
  expect_true(any(vapply(q$warnings, inherits, NA, "frmtmb_conv_code")))
})

test_that("an alias under a prior is not taken as flat", {
  # y ~ x + s(x): s(x).fx1 is x's exact alias, flat in the likelihood,
  # but a prior on the two curves the objective along it. With the
  # design's null directions taken as flat whatever the prior, stops short
  # along it under normal(0, 1e4) were called the maximum: code 8, 1.87
  # short, "so the estimates are at the maximum"; code 0 with no warning
  # (review of lane code8, fix-round check, dev/review-code8-fix-dnull.R).
  # The optimizer runs nlminb to the maximum and steps along the null
  set.seed(5)
  d <- data.frame(x = stats::runif(200, 0, 10))
  d$y <- 1 + 0.5 * d$x + sin(d$x) + stats::rnorm(200, 0, 0.5)
  pr <- set_prior("normal(0, 1e4)", class = "b", coef = "x") +
    set_prior("normal(0, 1e4)", class = "b", coef = "s(x).fx1")
  ref <- frm(bf(y ~ x + s(x)), data = d, prior = pr)
  v <- frmtmb:::se_design_null(ref)[, 1]
  step_off <- function(code) {
    function(par, fn, gr, lower, upper, control) {
      o <- stats::nlminb(par, fn, gr, control = list(iter.max = 1000,
                                                     eval.max = 2000))
      p <- o$par + 2e4 * v
      list(par = p, objective = fn(p), convergence = code,
           message = if (code == 1L) "false convergence (8)" else
             "relative convergence (4)")
    }
  }
  for (code in c(1L, 0L)) {
    r <- conv_classes(frm(bf(y ~ x + s(x)), data = d, prior = pr,
                          control = frmtmb_control(optimizer = step_off(code),
                                                   restarts = 0)))
    expect_gt(r$value$opt$objective - ref$opt$objective, 0.1)
    expect_false(any(vapply(r$messages, inherits, NA,
                            "frmtmb_converged_at_max")))
    wm <- vapply(r$warnings, conditionMessage, "")
    expect_false(any(grepl("so the estimates are at the maximum", wm,
                           fixed = TRUE)), info = paste("code", code))
    expect_true(any(vapply(r$warnings, inherits, NA,
                           c("frmtmb_conv_gradient", "frmtmb_conv_code"))),
                info = paste("code", code))
  }
})

conv_valley <- function(x0, seed) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(60, x0, x0 + 1))
  d$y <- 2 * exp((3 / x0) * d$x) + stats::rnorm(60, 0, 1)
  d
}
conv_valley_bf <- bf(y ~ A * exp(B * x), A ~ 1, B ~ 1, nl = TRUE)
# optim()'s BFGS cut at `maxit`, its code reported as it is or, with
# `fake8`, as nlminb's "false convergence (8)"
conv_bfgs <- function(maxit, fake8 = FALSE) {
  function(par, fn, gr, lower, upper, control) {
    o <- stats::optim(par, fn, gr, method = "BFGS",
                      control = list(maxit = maxit))
    list(par = o$par, objective = o$value,
         convergence = if (fake8) 1L else as.integer(o$convergence),
         message = if (fake8) "false convergence (8)" else "optim BFGS")
  }
}

test_that("a stop short on a curved valley is not called the maximum", {
  # y = A exp(B x) with x on [100, 101]: A and B trade off along a curved
  # valley. A first version of this round's convergence check took the
  # best of a few straight steps in place of the Newton gain, and called
  # these stops the maximum: optim()'s BFGS cut at 40 iterations and
  # reported as code 8, 0.262 short, got the classed message ("gains
  # 3.72e-06"), and BFGS's own code-1 stop 8.0 short got "so the
  # estimates are at the maximum" and no gradient warning (review of lane
  # code8, finding 1, dev/review-code8-valley8.R and -valley-opt.R)
  tol <- frmtmb_control()$grad_tol
  d1 <- conv_valley(100, 1)
  ref1 <- frm(conv_valley_bf, data = d1)
  r <- conv_classes(frm(conv_valley_bf, data = d1,
                        control = frmtmb_control(
                          optimizer = conv_bfgs(40, fake8 = TRUE),
                          restarts = 0)))
  expect_gt(r$value$opt$objective - ref1$opt$objective, 10 * tol)
  expect_false(any(vapply(r$messages, inherits, NA,
                          "frmtmb_converged_at_max")))
  wm <- vapply(r$warnings, conditionMessage, "")
  expect_true(any(vapply(r$warnings, inherits, NA, "frmtmb_conv_code")))
  expect_false(any(grepl("so the estimates are at the maximum", wm,
                         fixed = TRUE)))
  d2 <- conv_valley(100, 2)
  ref2 <- frm(conv_valley_bf, data = d2)
  r <- conv_classes(frm(conv_valley_bf, data = d2,
                        control = frmtmb_control(optimizer = conv_bfgs(5000),
                                                 restarts = 0)))
  # (BFGS's own stop: code 1 here; a platform where it reaches the
  # maximum has nothing to report)
  if (r$value$opt$objective - ref2$opt$objective > 10 * tol) {
    wm <- vapply(r$warnings, conditionMessage, "")
    expect_false(any(grepl("so the estimates are at the maximum", wm,
                           fixed = TRUE)))
    expect_true(any(vapply(r$warnings, inherits, NA,
                           c("frmtmb_conv_gradient", "frmtmb_conv_code"))))
  }
})

test_that("a stop at a correlation's limit is measured as the maximum", {
  # brms's example fit 1 (student, mo(), s(), arma, sigma ~ Trt) on its
  # 40 rows: the correlation of (1 + Trt | visit) runs to -1 (its
  # parameter at -723 to -4130, 4 visits), the smooth's sd to 0, and
  # s(Age)'s fixed column is Age's exact alias. Up to 0.80.0 the
  # correlation's row left the curvature at the stop unusable, so the
  # fit was called short and kept sdreport()'s inverse, with negative
  # variances for every parameter: fitted(dpar = "sigma") had NaN errors
  # (ported row brmsfit-methods:326). The boundary is held, as the
  # standard-error check holds it, and the rest is judged as a
  # converged fit is (dev/code8-findings.md, Part B)
  skip_if_not_installed("brms")
  ex <- get("brmsfit_example1", envir = asNamespace("brms"))
  d <- as.data.frame(ex$data)
  attr(d, "terms") <- NULL
  r <- conv_classes(frm(bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume +
                             offset(Age) + (1 + Trt | visit) +
                             arma(visit, patient), sigma ~ Trt),
                          family = student(), data = d))
  fit <- r$value
  wm <- vapply(r$warnings, conditionMessage, "")
  expect_false(any(grepl("Large maximum absolute gradient", wm)))
  expect_false(any(grepl("short of the maximum|curvature at the stop", wm)))
  bnd <- Filter(function(x) inherits(x, "frmtmb_boundary_fit"), r$messages)
  expect_length(bnd, 1L)
  expect_match(conditionMessage(bnd[[1L]]), "cor_visit__Intercept__Trt1",
               fixed = TRUE)
  lost <- sdr_of(fit)$se_lost
  # theta_3, the block's correlation parameter
  expect_identical(unname(lost[["theta_3"]]), "boundary")
  se <- sqrt(diag(allow_warnings(vcov(fit, full = TRUE),
                                 "Standard errors are not available")))
  nm <- frmtmb:::outer_par_names(fit)
  expect_true(all(is.finite(se[match(c("sigma_(Intercept)", "sigma_Trt1"),
                                     nm)])))
  # the alias is named, and a coefficient off it keeps its SE
  expect_true(all(c("Age", "s(Age).fx1") %in% names(lost)))
  expect_true(is.finite(se[[match("volume", nm)]]))
  fi <- allow_warnings(fitted(fit, dpar = "sigma"),
                       "Standard errors are not available")
  expect_true(all(is.finite(fi)) && all(fi > 0))
})

test_that("a converged fit without a convergence warning is unchanged", {
  # the guard's absent case: the clause belongs to the convergence
  # warning only
  d <- conv_alias_data()
  fit <- frm(y ~ x, data = d)
  r <- conv_conds(vcov(fit))
  expect_length(r$warnings, 0L)
  expect_equal(fit$opt$convergence, 0L)
})

test_that("an sd held near zero short of its maximum is not at the maximum", {
  # The sd's row is empty and its gradient is 5.6e-10, so a Newton step
  # cannot see what a larger sd gains. The headroom asks the likelihood
  # (se_sd_gain_up()); without that, this stop 92 below the maximum was
  # called the maximum (review of lane illse, round 2, R2-1). The fit
  # starts at lm()'s maximum given sd = 0 and takes no step, so the sd
  # is the only direction that is short whatever start frm() would pick;
  # from a start where another coefficient is short too, the stop is
  # called short without R2-1 (the 0.71.0 release review)
  set.seed(1)
  n <- 300
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(rep(1:20, length.out = n)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.5)[d$g] +
    stats::rnorm(n, 0, 0.4)
  ls0 <- stats::lm(y ~ x, data = d)
  start <- list(beta = unname(stats::coef(ls0)),
                betad = log(sqrt(mean(stats::residuals(ls0)^2))),
                theta = -15)
  r <- conv_conds(frm(y ~ x + (1 | g), data = d, start = start,
                      control = frmtmb_control(optCtrl = list(
                        iter.max = 0, eval.max = 1000))))
  expect_false(r$value$opt$convergence == 0L)
  w <- grep(conv_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(w, 1L)
  expect_match(w, "short of the maximum", fixed = TRUE)
  expect_false(any(grepl("so the estimates are at the maximum",
                         r$warnings, fixed = TRUE)))
})
