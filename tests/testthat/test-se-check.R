# Standard errors the outer Hessian cannot give (R/se-check.R): the
# repairs that recover them, the fit-time warning that names the rest,
# and summary()'s reason line. Lane nanse, dev/nanse-findings.md.

# brms_monotonic's own data code, as dev/nanse-mo-sweep.R builds it
se_mo_data <- function(seed) {
  set.seed(seed)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}

se_lost_phrase <- "Standard errors are not available"

test_that("a mo() simplex on its boundary keeps a finite Hessian", {
  # seed 7: the interaction's simplex has a weight at 0. On the softmax
  # coordinates of 0.68.1 that weight sat at a coordinate near -46, the
  # outer Hessian's reciprocal condition number was 3.7e-25, solve()
  # refused it, and the SE check held the coordinate fixed to get the
  # rest. The simplex chart of mo_simplex() reaches the face at a
  # finite coordinate with a finite curvature (lane optima), so the
  # plain inverse exists and is what the fit reports.
  d <- se_mo_data(7)
  w <- character()
  fit <- withCallingHandlers(frm(ls ~ mo(income) * age, data = d),
                             warning = function(x) {
                               w <<- c(w, conditionMessage(x))
                               invokeRestart("muffleWarning")
                             })
  expect_identical(fit$opt$convergence, 0L)
  expect_length(w, 0L)
  H <- fit$obj$he(fit$opt$par)
  V <- tryCatch(solve(H), error = function(e) NULL)
  expect_false(is.null(V))
  expect_true(all(is.finite(V)) && all(diag(V) > 0))
  # the fixture still has its weight at the face
  s <- frmtmb:::mo_simplex(fit$estimates$zeta2)
  expect_lt(min(s) / max(s), 1e-8)
  se <- fixef(fit)[, "Est.Error"]
  expect_equal(unname(se), sqrt(diag(V))[seq_len(4L)], tolerance = 1e-4)
  # the vignette's next call printed "need finite 'ylim' values"
  ce <- conditional_effects(fit, "income:age")
  expect_true(all(is.finite(ce[[1L]]$se__)))
})

test_that("a lost simplex SE warns at fit time and keeps the rest", {
  # helper-mo-flat.R: the simplex coordinate does not enter the
  # likelihood, so its row of the Hessian is exactly zero and no
  # inverse of the whole matrix exists
  w <- character()
  fit <- withCallingHandlers(frm(ls ~ mo(inc) + age, data = mo_flat_data()),
                             warning = function(x) {
                               w <<- c(w, conditionMessage(x))
                               invokeRestart("muffleWarning")
                             })
  expect_identical(fit$opt$convergence, 0L)
  expect_length(w, 1L)
  expect_match(w, se_lost_phrase, fixed = TRUE)
  expect_match(w, "1 of 5 parameters", fixed = TRUE)
  expect_match(w, "zeta1_1 (simplex of moinc)", fixed = TRUE)
  expect_match(w, "flat", fixed = TRUE)
  expect_true(all(is.finite(fixef(fit)[, "Est.Error"])))
  ci <- confint(fit)
  expect_true(is.nan(ci["zeta1_1", "lwr"]))
  expect_true(all(is.finite(ci[rownames(ci) != "zeta1_1", "lwr"])))
  out <- capture.output(print(summary(fit)))
  expect_true(any(grepl("without a standard error", out, fixed = TRUE)))
  expect_true(any(grepl("zeta1_1", out, fixed = TRUE)))
  # said once: vcov() does not repeat it
  expect_silent(vcov(fit))
  # a prediction does not read the lost weight's direction, so it
  # keeps its standard error
  ce <- conditional_effects(fit, "inc:age")
  expect_true(all(is.finite(ce[[1L]]$se__)))
})

test_that("a covariate spanning 1e5 in a nonlinear body keeps its SEs", {
  # optimHess() steps 1e-3 in b = -2e-5, which moves b * x by 100: the
  # Hessian's b diagonal came out 1.3e96 and every SE was NaN, in
  # silence. The exact Hessian gives the SEs of the rescaled control.
  set.seed(41)
  x <- runif(300, 0, 1e5)
  d <- data.frame(x = x, xs = x / 1e5,
                  y = 3 * exp(-2e-5 * x) + rnorm(300, 0, 0.05))
  f1 <- frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), data = d,
            start = list(beta = c(3, -2e-5)))
  f2 <- frm(bf(y ~ a * exp(b * xs), a ~ 1, b ~ 1, nl = TRUE), data = d,
            start = list(beta = c(3, -2)))
  s1 <- fixef(f1)[, "Est.Error"]
  s2 <- fixef(f2)[, "Est.Error"]
  expect_true(all(is.finite(s1)))
  expect_equal(unname(s1), unname(s2 * c(1, 1e-5)), tolerance = 1e-4)
})

test_that("an unidentified ridge names its parameters and keeps sigma", {
  # a + b + log(c0) with a, b ~ 1 + x: only a + b + log(c0) and
  # a_x + b_x are identified. With c0 at 5.9e-5 optimHess() steps into
  # c0 < 0, its Hessian row is NaN, and lane fixes' flat-direction
  # check returns NULL there; the fit used to say nothing.
  set.seed(45)
  d <- data.frame(x = rnorm(200))
  d$y <- 1 + 0.5 * d$x + log(5e-5) + rnorm(200, 0, 0.3)
  w <- character()
  f <- withCallingHandlers(
    frm(bf(y ~ a + b + log(c0), a ~ 1 + x, b ~ 1 + x, c0 ~ 1, nl = TRUE),
        data = d, start = list(beta = c(0.5, 0.25, 0.5, 0.25, 5e-5))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  hit <- grep(se_lost_phrase, w, fixed = TRUE, value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, "5 of 6 parameters", fixed = TRUE)
  for (p in c("a_(Intercept)", "a_x", "b_(Intercept)", "b_x",
              "c0_(Intercept)")) {
    expect_match(hit, p, fixed = TRUE)
  }
  # sigma is identified and gets the standard error of the same mean
  # written without the ridge
  f0 <- frm(y ~ x, data = d)
  se_sig <- function(fit) {
    vc <- vcov(fit, full = TRUE)
    sqrt(vc["sigma_(Intercept)", "sigma_(Intercept)"])
  }
  expect_equal(se_sig(f), se_sig(f0), tolerance = 1e-4)
  ci <- confint(f)
  expect_true(all(is.nan(ci[c("a_(Intercept)", "b_x"), "lwr"])))
})

test_that("a bound-held parameter loses its SE and the others keep theirs", {
  # test-backlog "The covariance machinery is not bound-aware":
  # gradcheck-01's construction, seed 101. The full Hessian is
  # indefinite at the constrained optimum and every SE was NaN.
  set.seed(101)
  dd <- data.frame(x = rnorm(300))
  dd$y <- rnorm(300, 1 + 2 * dd$x, 1)
  w <- character()
  f <- withCallingHandlers(
    frm(bf(y ~ x), family = gaussian(), data = dd,
        prior = set_prior("", class = "b", ub = 0.1)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  hit <- grep(se_lost_phrase, w, fixed = TRUE, value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, "x: a bound holds", fixed = TRUE)
  V <- vcov(f, full = TRUE)
  nm <- colnames(V)
  free <- nm != "x"
  expect_true(is.nan(V["x", "x"]))
  H <- f$obj$he(f$opt$par)
  ref <- solve(H[free, free])
  expect_equal(unname(V[free, free]), unname(ref), tolerance = 1e-6)
})

test_that("a healthy fit is untouched and never warns", {
  # the guard's absent case: sdreport()'s own inverse is usable, so the
  # covariance is exactly sdreport()'s
  d <- se_mo_data(7)
  expect_silent(fit <- frm(ls ~ mo(income), data = d))
  V <- vcov(fit, full = TRUE)
  ref <- RTMB::sdreport(fit$obj)$cov.fixed
  expect_identical(unname(V), unname(ref))
})

test_that("check_se = 'ignore' and 'stop' do what they say", {
  d <- mo_flat_data()
  expect_silent(frm(ls ~ mo(inc) + age, data = d,
                    control = frmtmb_control(check_se = "ignore")))
  expect_error(frm(ls ~ mo(inc) + age, data = d,
                   control = frmtmb_control(check_se = "stop")),
               se_lost_phrase, fixed = TRUE)
})

test_that("se = TRUE says it once, and a fit that stopped short not at all", {
  d <- mo_flat_data()
  w <- character()
  withCallingHandlers(frm(ls ~ mo(inc) + age, data = d, se = TRUE),
                      warning = function(x) {
                        w <<- c(w, conditionMessage(x))
                        invokeRestart("muffleWarning")
                      })
  expect_length(grep(se_lost_phrase, w, fixed = TRUE), 1L)
  expect_length(w, 1L)
  # an optimizer stopped by its iteration cap: the convergence warning
  # is the verdict, and the curvature there is not worth a second one
  w <- character()
  short <- withCallingHandlers(
    frm(ls ~ mo(inc) + age, data = d,
        control = frmtmb_control(optCtrl = list(iter.max = 3,
                                                eval.max = 5),
                                 restarts = 0)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("did not report convergence", w, fixed = TRUE)))
  expect_length(grep(se_lost_phrase, w, fixed = TRUE), 0L)
  # and its covariance is sdreport()'s own, unrepaired
  V <- suppressWarnings(vcov(short, full = TRUE))
  ref <- suppressWarnings(RTMB::sdreport(short$obj)$cov.fixed)
  expect_identical(unname(V), unname(ref))
})

# ---- punch round 1 (dev/reviews/2026-10-06-nanse.md) ------------------

# the review's spread ridge: y ~ a + b, a ~ 0 + f, b ~ 1. Only a_k + b
# is identified; the null vector loads sqrt(1 / (2k)) on each a_k
se_spread <- function(k, m = 5, seed = 1) {
  set.seed(seed)
  f <- factor(rep(seq_len(k), each = m))
  data.frame(f = f, y = rnorm(k)[f] + rnorm(k * m, 0, 0.5))
}

se_capture <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

test_that("a flat direction shared by many coefficients takes all of them", {
  # B1: at k = 60 each a_k loads 0.091 on the flat direction, under the
  # old loading threshold of 0.1, and kept an SE of 0.197 with z = 7.1
  r <- se_capture(frm(bf(y ~ a + b, a ~ 0 + f, b ~ 1, nl = TRUE),
                      data = se_spread(60)))
  # on the merged tree lane fixes' flat-direction warning speaks first and
  # explains the lost SEs; it must name all 60 a_f* coefficients
  hit <- grep("are not identified: at the optimum", r$warnings,
              fixed = TRUE, value = TRUE)
  expect_length(hit, 1L)
  expect_length(grep(se_lost_phrase, r$warnings, fixed = TRUE), 0L)
  expect_length(regmatches(hit, gregexpr("a_f[0-9]+", hit))[[1]], 60L)
  se <- fixef(r$value)[, "Est.Error"]
  expect_true(all(is.nan(se)))
  ci <- suppressWarnings(confint(r$value))
  expect_identical(sum(is.finite(ci[, "lwr"])), 1L)   # sigma alone
  # two crossed factors: every coefficient unidentified
  set.seed(3)
  k <- 20
  d <- expand.grid(f = factor(seq_len(k)), g = factor(seq_len(k)))
  d <- d[sample(nrow(d), 6 * k), ]
  d$y <- rnorm(k)[d$f] + rnorm(k)[d$g] + rnorm(nrow(d), 0, 0.5)
  r <- se_capture(frm(bf(y ~ a + b, a ~ 0 + f, b ~ 0 + g, nl = TRUE),
                      data = d))
  expect_true(all(is.nan(fixef(r$value)[, "Est.Error"])))
})

test_that("a small loading on a downward direction keeps its SE", {
  # B1's other half and m1: frmtmb.sample's gr(g, by = f) fixture. A
  # slope SD at exp(-6.6) and its correlation load 0.70 each on an
  # eigenvalue of -1.18, x loads 0.10, and a step along it gains 1e-4 in
  # log-likelihood. x keeps an SE (lme4 0.1009 on the same model);
  # the two named parameters are flat, not "not a maximum".
  set.seed(11)
  dd <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, 10)))
  dd$f <- factor(ifelse(as.integer(dd$g) <= 8, "a", "b"))
  u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
  dd$y <- stats::rnorm(160, 1 + 0.5 * dd$x + u[dd$g, 1] +
                         u[dd$g, 2] * dd$x, 1)
  r <- se_capture(frm(bf(y ~ x + (1 + x | gr(g, by = f))),
                      family = gaussian(), data = dd))
  hit <- grep(se_lost_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, "theta_2, theta_3: the likelihood is flat", fixed = TRUE)
  expect_false(grepl("not a maximum", hit, fixed = TRUE))
  se_x <- fixef(r$value)["x", "Est.Error"]
  # the plain sdreport() inverse, which keeps x finite on this fit
  raw <- suppressWarnings(RTMB::sdreport(r$value$obj)$cov.fixed)
  expect_equal(se_x, sqrt(raw[2, 2]), tolerance = 0.05)
})

test_that("a prediction along a lost direction gets NaN and one warning", {
  # B2: round 0 of the lane gave fitted(dpar = "a") an Est.Error of
  # 0.1818 from the pseudo-inverse, silently
  fit <- suppressWarnings(frm(bf(y ~ a + b, a ~ 0 + f, b ~ 1, nl = TRUE),
                              data = se_spread(10)))
  nd <- se_spread(10)[c(1, 6), ]
  for (dp in c("a", "b")) {
    r <- se_capture(fitted(fit, newdata = nd, dpar = dp))
    expect_true(all(is.nan(r$value[, "Est.Error"])))
    expect_true(any(grepl("predictions move along a direction",
                          r$warnings, fixed = TRUE)))
  }
  r <- se_capture(frm_linpred(fit, newdata = nd, dpar = "a",
                              se.fit = TRUE))
  expect_true(all(is.nan(r$value$se.fit)))
  expect_length(r$warnings, 1L)
  # the absent case: a lost simplex coordinate the prediction does not
  # read (its gradient is exactly zero) keeps every band
  fm <- suppressWarnings(frm(ls ~ mo(inc) + age, data = mo_flat_data()))
  r <- se_capture(conditional_effects(fm, "inc:age"))
  expect_true(all(is.finite(r$value[[1L]]$se__)))
  expect_false(any(grepl("predictions move", r$warnings, fixed = TRUE)))
})

test_that("hypothesis() and emmeans() treat a lost direction the same way", {
  # B2: hypothesis() read the masked covariance, so a_f1 + b, which the
  # data determine, got NaN; emmeans() gave NaN without the warning
  dd <- se_spread(10)
  fit <- suppressWarnings(frm(bf(y ~ a + b, a ~ 0 + f, b ~ 1, nl = TRUE),
                              data = dd))
  # a_f1 + b is the mean of cell 1, whose ML standard error is the ML
  # residual sd over sqrt(5)
  s_ml <- sqrt(mean((dd$y - ave(dd$y, dd$f))^2))
  r <- se_capture(hypothesis(fit, "a_f1 + b_Intercept = 0"))
  expect_equal(r$value$hypothesis$Est.Error, s_ml / sqrt(5),
               tolerance = 1e-5)
  expect_length(r$warnings, 0L)
  r <- se_capture(hypothesis(fit, c("a_f1 = 0", "a_f1 + b_Intercept = 0")))
  expect_true(is.nan(r$value$hypothesis$Est.Error[1]))
  expect_true(is.finite(r$value$hypothesis$Est.Error[2]))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "1 of 2 hypotheses move along", fixed = TRUE)
  skip_if_not_installed("emmeans")
  r <- se_capture(summary(emmeans::emmeans(fit, ~ f, dpar = "a")))
  expect_true(all(is.nan(r$value$SE)))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "10 of 10 predictions move along", fixed = TRUE)
  r <- se_capture(summary(emmeans::emmeans(fit, ~ f)))
  expect_equal(r$value$SE, rep(s_ml / sqrt(5), 10), tolerance = 1e-5)
  expect_length(r$warnings, 0L)
})

test_that("an earlier warning explains only the parameters it names", {
  # m3: a one-level (1 | g1) warns about its own variance; x, held by a
  # bound, still gets the SE warning
  set.seed(101)
  d3 <- data.frame(x = rnorm(300), g1 = factor("only"))
  d3$y <- rnorm(300, 1 + 2 * d3$x, 1)
  r <- se_capture(frm(bf(y ~ x + (1 | g1)), family = gaussian(), data = d3,
                      prior = set_prior("", class = "b", ub = 0.1)))
  hit <- grep(se_lost_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, "1 of 4 parameters. x: a bound holds it", fixed = TRUE)
  expect_true(any(grepl("single level", r$warnings, fixed = TRUE)))
  # and vcov() does not follow with an older, different reason
  expect_silent(vcov(r$value))
  # check_olre = "ignore" silences the OLRE warning, so it explains
  # nothing, and the SE warning names the two parameters
  set.seed(5)
  d2 <- data.frame(id = factor(1:80), x = rnorm(80))
  d2$y <- 1 + 0.5 * d2$x + rnorm(80)
  r <- se_capture(frm(bf(y ~ x + (1 | id)), family = gaussian(), data = d2,
                      control = frmtmb_control(check_olre = "ignore")))
  hit <- grep(se_lost_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, "sigma_(Intercept), theta_1", fixed = TRUE)
  # with the OLRE warning on, it explains both, and nothing follows
  r <- se_capture(frm(bf(y ~ x + (1 | id)), family = gaussian(), data = d2))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "confounded with the residual sd", fixed = TRUE)
  expect_silent(vcov(r$value))
})

test_that("REML does not promise the coefficients' standard errors", {
  # m5: under REML the coefficients come from the joint precision, which
  # the repair does not reach
  set.seed(5)
  d2 <- data.frame(id = factor(1:80), x = rnorm(80))
  d2$y <- 1 + 0.5 * d2$x + rnorm(80)
  r <- se_capture(frm(bf(y ~ x + (1 | id)), family = gaussian(), data = d2,
                      REML = TRUE,
                      control = frmtmb_control(check_olre = "ignore")))
  hit <- grep(se_lost_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, "joint precision", fixed = TRUE)
  expect_false(grepl("other standard errors are kept", hit, fixed = TRUE))
})

test_that("separated data are named as separation", {
  # m7: test-predfix's design, seed 514. The estimates run off toward
  # infinity, so where nlminb stops is a matter of rounding: at the
  # default budget of 1000 it stopped at code 0 on Windows and hit the
  # limit (code 9) on Ubuntu, whose convergence warning then explains
  # the fit and the SE check stays silent. At 4000 it stops at code 0
  # on both, after 1939 objective and gradient calls here and 2090 with
  # OpenBLAS 0.3.26 (dev/cifix-findings.md).
  set.seed(514)
  d <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
  d$yb <- as.integer(d$z > 0)
  ctl <- frmtmb_control(optCtrl = list(eval.max = 4000, iter.max = 4000))
  r <- se_capture(frm(yb ~ z + x, family = bernoulli(), data = d,
                      control = ctl))
  expect_identical(r$value$opt$convergence, 0L)
  hit <- grep(se_lost_phrase, r$warnings, fixed = TRUE, value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, "the data separate the outcomes", fixed = TRUE)
})

test_that("the fit-time Hessian is sdreport()'s", {
  # m2: fd_hessian() is optimHess() step for step, so sdreport() is
  # handed the matrix it would build, names included
  fit <- frm(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy)
  obj <- fit$obj
  ns <- asNamespace("frmtmb")
  st <- ns$obj_state_save(obj)
  at <- ns$sdr_outer_point(obj)
  a <- stats::optimHess(at, obj$fn, obj$gr)
  ns$obj_state_restore(obj, st)
  b <- ns$fd_hessian(at, obj$gr)$H
  ns$obj_state_restore(obj, st)
  expect_identical(a, b)
})

test_that("whether the check waits is decided by counted work", {
  # RB1: a wall-clock budget deferred the same call 1 time in 20 idle
  # and 3 in 20 under load. The decision reads the parameter count and
  # the optimizer's counted evaluations, nothing else.
  ns <- asNamespace("frmtmb")
  at_fit <- function(np, ev) {
    ns$se_check_at_fit(list(opt = list(par = numeric(np), evals = ev)))
  }
  expect_true(at_fit(10, NULL))
  expect_false(at_fit(50, NULL))
  expect_false(at_fit(50, 399))
  expect_true(at_fit(50, 400))
  fits <- lapply(1:2, function(i) {
    frm(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy)
  })
  expect_gt(fits[[1]]$opt$evals, 0)
  expect_identical(fits[[1]]$opt$evals, fits[[2]]$opt$evals)
})

test_that("a waiting check reaches the caller, even through fixef()", {
  # RB1: fixef() reads vcov() under suppressWarnings(); the report was
  # raised there, muffled, and marked as given
  set.seed(5)
  d2 <- data.frame(id = factor(1:80), x = rnorm(80))
  d2$y <- 1 + 0.5 * d2$x + rnorm(80)
  waiting <- function() {
    f <- frm(bf(y ~ x + (1 | id)), family = gaussian(), data = d2,
             control = frmtmb_control(check_olre = "ignore",
                                      check_se = "ignore"))
    f$cache$se_deferred <- "warning"
    f
  }
  f2 <- waiting()
  r <- se_capture(fixef(f2))
  expect_length(grep(se_lost_phrase, r$warnings, fixed = TRUE), 1L)
  expect_silent(vcov(f2))
  # under a suppressWarnings() that frmtmb wrote it stays pending
  f3 <- waiting()
  quiet <- function(x) suppressWarnings(vcov(x))
  environment(quiet) <- asNamespace("frmtmb")
  quiet(f3)
  r <- se_capture(vcov(f3))
  expect_length(grep(se_lost_phrase, r$warnings, fixed = TRUE), 1L)
  expect_silent(vcov(f3))
  # a user's own suppressWarnings() is the user declining it
  f4 <- waiting()
  suppressWarnings(vcov(f4))
  expect_silent(vcov(f4))
})

test_that("a lost group sd does not take the other variance components", {
  # RB2: VarCorr() multiplied the shown covariance, whose lost rows are
  # NaN, so g1's sd and the residual sd lost their errors too
  set.seed(1)
  n <- 60
  d <- data.frame(x = rnorm(n), g1 = factor(rep(1:12, 5)), g2 = gl(6, 10))
  d$y <- 1 + 0.5 * d$x + rnorm(12, 0, 0.8)[d$g1] + rnorm(n)
  fit <- suppressWarnings(frm(bf(y ~ x + (1 | g1) + (1 + x | g2)),
                              data = d))
  expect_setequal(names(sdr_of(fit)$se_lost),
                  c("theta_2", "theta_3", "theta_4"))
  r <- se_capture(VarCorr(fit))
  vc <- r$value
  h <- suppressWarnings(hypothesis(fit, "sd_g1__Intercept = 0",
                                   class = NULL))
  expect_equal(vc$g1$sd[1, "Est.Error"], h$hypothesis$Est.Error,
               tolerance = 1e-6)
  expect_true(is.finite(vc$residual__$sd[1, "Est.Error"]))
  expect_true(all(is.nan(vc$g2$sd[, "Est.Error"])))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "VarCorr() entries move along", fixed = TRUE)
  s <- summary(fit)
  expect_true(is.finite(s$random$g1[1, "Est.Error"]))
  expect_true(is.finite(s$spec_pars[1, "Est.Error"]))
  r <- se_capture(conditional_effects(fit, "x", re_formula = NULL))
  expect_true(all(is.finite(r$value[[1L]]$se__)))
})

test_that("a lost simplex coordinate does not take the other predictions", {
  # m1 follow-up: on mo() seed 12 a concave direction of zeta2_2 loaded
  # 0.022 on zeta1_2, which kept its SE; the whole vector in the null
  # basis took every conditional_effects() band, and the plot stopped
  # on "need finite 'ylim' values". That fit stood on a softmax plateau
  # short of the maximum, which the simplex chart of lane optima no
  # longer gives (dev/optima-findings.md), so the fixture is now two
  # exactly flat coordinates, one per term (helper-mo-flat.R).
  fm <- suppressWarnings(frm(ls ~ mo(inc) * age, data = mo_flat_data()))
  expect_identical(unname(sdr_of(fm)$se_lost), c("flat", "flat"))
  r <- se_capture(conditional_effects(fm, "age"))
  expect_true(all(is.finite(r$value[[1L]]$se__)))
  expect_length(r$warnings, 0L)
  # a display along the lost direction draws its line without a band
  ce <- suppressWarnings(conditional_effects(fm, "inc:age"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(ce, ask = FALSE))
})

test_that("an identified remainder keeps its SEs after a flat row goes", {
  # RB4: (1 + x | g2) with no g2 variation empties theta_3's row; the
  # rest is positive definite, but its smallest eigenvalue (0.022 on
  # seed 11, the intercept against the 39 contrasts) sat below a noise
  # threshold summed over the whole 44 x 44 block, and the intercept and
  # f2..f40 lost SEs that lme4 reports
  skip_if_not_installed("lme4")
  for (s in c(11L, 17L)) {
    set.seed(s)
    n <- 480
    d <- data.frame(x = rnorm(n), f = factor(sample(1:40, n, TRUE)),
                    g2 = factor(rep(1:6, length.out = n)))
    d$y <- 1 + 0.5 * d$x + rnorm(40, 0, 0.5)[d$f] + rnorm(n)
    fit <- suppressWarnings(frm(y ~ x + f + (1 + x | g2), data = d))
    # 45 outer parameters: the check waited, and reports here
    lost <- names(suppressWarnings(sdr_of(fit))$se_lost)
    expect_true(all(startsWith(lost, "theta_")), info = paste("seed", s))
    m <- suppressMessages(suppressWarnings(
      lme4::lmer(y ~ x + f + (1 + x | g2), data = d, REML = FALSE)))
    ref <- sqrt(diag(as.matrix(stats::vcov(m))))[c("(Intercept)", "f2")]
    se <- suppressWarnings(sqrt(diag(vcov(fit))))[c("Intercept", "f2")]
    # base sdreport() and lme4 differ by 0.1 to 0.4 percent here (one
    # sits on the boundary, the other at a tiny sd); 2 percent of the SE
    # itself separates that from a lost or a meaningless value
    expect_equal(unname(se), unname(ref), tolerance = 0.02,
                 info = paste("seed", s))
  }
})

test_that("a fit with random effects and a lost sd gets a covariance", {
  # Both smoothing sds run to zero and lose their standard errors, so
  # the outer Hessian is indefinite and so is the joint precision built
  # on it. 0.68.0 inverted that precision as it was: grid rows got
  # coefficient variances as low as -0.0026 here, and
  # frm_linpred(se.fit = TRUE) reported the kriging variance alone
  # (dev/cifix-findings.md, scan seed 21, n = 100). The repaired
  # covariance propagates only the outer parameters that keep a
  # standard error, so no row can fall below its variance with every
  # outer parameter held at its estimate: that is the law of total
  # variance, and it holds on any fit, lost parameters or not.
  set.seed(21)
  n <- 100
  d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                  fac = factor(rep(c("A", "B"), length.out = n)))
  d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
    stats::rnorm(n, 0, 0.3)
  fit <- suppressWarnings(
    frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)), data = d))
  gx <- d$x[-1] - diff(d$x) / 2
  nd <- data.frame(x = rep(gx, 2),
                   fac = factor(rep(c("A", "B"), each = length(gx))))
  lb <- frm_lp_basis(fit, newdata = nd)
  A <- as.matrix(lb$A)
  q <- rowSums((A %*% lb$V) * A)
  Q <- joint_precision(fit)
  r <- fit$obj$env$random
  pr <- match(lb$coef_pos, r)
  inr <- !is.na(pr)
  W <- as.matrix(Matrix::solve(Q[r, r]))[pr[inr], pr[inr]]
  q_in <- rowSums((A[, inr] %*% W) * A[, inr])
  expect_true(all(q_in > 0))
  expect_gt(min(q / q_in), 1 - 1e-8)
  # and frm_linpred() reports that variance, not the kriging part alone
  se <- frm_linpred(fit, newdata = nd, se.fit = TRUE)$se.fit
  expect_equal(se^2, q + lb$extra_var, tolerance = 1e-10)
  # The rest depends on this platform's fit having lost the sds, which
  # it does on Windows with the reference BLAS and with OpenBLAS 0.3.26;
  # a run where it did not has asserted the bound above and stops here.
  lost <- sdr_of(fit)$se_lost
  skip_if(!length(lost), "this fit kept every standard error here")
  # frm_joint_cov() shows a lost parameter as vcov() does, NaN, and the
  # covariance of everything else is finite
  jc <- frm_joint_cov(fit)
  bad <- jc$lost_pos
  expect_length(bad, length(lost))
  expect_true(all(jc$names[bad] == "theta"))
  expect_true(all(is.nan(jc$V[bad, ])))
  expect_true(all(is.finite(jc$V[-bad, -bad])))
})
