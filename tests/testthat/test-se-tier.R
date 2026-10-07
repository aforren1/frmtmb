# Standard errors frmtmb used to report from noise, and their names:
# the standard-error check's tiers ask tier 3 first (R/se-check.R,
# cov_from_hessian()), a variance component at zero is a boundary fit,
# separation is named whatever the optimizer's code, and ranef()'s
# conditional SDs read the repaired joint covariance. Lane setier,
# dev/setier-findings.md.

tier_conds <- function(expr) {
  w <- character()
  m <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m[[length(m) + 1L]] <<- x
    invokeRestart("muffleMessage")
  })
  list(value = v, warnings = w, messages = m)
}

tier_sleep_nested <- function() {
  d <- lme4::sleepstudy
  d$a <- factor(d$Days %% 3)
  d
}

test_that("a variance at zero loses its SE whatever the sign of noise", {
  skip_if_not_installed("lme4")
  # Subject:a's sd sits at exp(-29): its row of the finite-difference
  # Hessian is +7.1e-12 with the reference BLAS and -7.1e-12 with
  # OpenBLAS, so solve() gave it a standard error of 375,150 on one
  # platform and NaN, with a warning, on the other (dev/cifix-diagnest.R)
  r <- tier_conds(frm(Reaction ~ Days + (1 | Subject/a),
                      data = tier_sleep_nested()))
  fit <- r$value
  expect_length(r$warnings, 0L)
  bm <- Filter(function(x) inherits(x, "frmtmb_boundary_fit"), r$messages)
  expect_length(bm, 1L)
  expect_match(conditionMessage(bm[[1L]]),
               "Boundary (singular) fit: sd_Subject:a__Intercept is at zero",
               fixed = TRUE)
  lost <- sdr_of(fit)$se_lost
  expect_identical(names(lost), "theta_1")
  expect_identical(unname(lost), "boundary")
  # the same verdict for either sign of the noise row
  ns <- asNamespace("frmtmb")
  h <- fit$cache$hessian_fixed
  j <- match("theta_1", ns$outer_par_names(fit))
  verdict <- function(sgn) {
    H <- h$H
    H[j, ] <- H[, j] <- sgn * abs(H[j, ])
    an <- ns$cov_from_hessian(fit, H, h$E)
    list(lost = an$lost, V = an$V[-j, -j])
  }
  pos <- verdict(1)
  neg <- verdict(-1)
  expect_identical(pos$lost, c(theta_1 = "flat"))
  expect_identical(neg$lost, pos$lost)
  # the others keep the covariance with that row held out
  ref <- solve(h$H[-j, -j])
  expect_equal(pos$V, ref, tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(neg$V, ref, tolerance = 1e-8, ignore_attr = TRUE)
  se <- sqrt(diag(vcov(fit, full = TRUE)))
  expect_true(is.nan(se[[j]]))
  expect_true(all(is.finite(se[-j])))
  # VarCorr() gives the sd no error and says nothing more about it
  v <- tier_conds(VarCorr(fit))
  expect_length(v$warnings, 0L)
  expect_true(is.nan(v$value[["Subject:a"]]$sd[1, "Est.Error"]))
  expect_true(is.finite(v$value$Subject$sd[1, "Est.Error"]))
  # summary() says why, and diagnose() lists the component
  expect_true(any(grepl("a variance at zero", summary(fit)$se_lost)))
  expect_false(is.null(diagnose(fit, quiet = TRUE)$singular))
})

test_that("check_se decides the boundary message too", {
  skip_if_not_installed("lme4")
  d <- tier_sleep_nested()
  r <- tier_conds(frm(Reaction ~ Days + (1 | Subject/a), data = d,
                      control = frmtmb_control(check_se = "ignore")))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 0L)
  expect_error(frm(Reaction ~ Days + (1 | Subject/a), data = d,
                   control = frmtmb_control(check_se = "stop")),
               "Boundary (singular) fit", fixed = TRUE)
})

test_that("a healthy fit with random effects is untouched and silent", {
  skip_if_not_installed("lme4")
  # the guard's absent case: no empty row and no flat direction, so the
  # covariance is sdreport()'s own, bit for bit, and nothing is said
  r <- tier_conds(frm(Reaction ~ Days + (Days | Subject),
                      data = lme4::sleepstudy))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 0L)
  expect_length(sdr_of(r$value)$se_lost, 0L)
  ref <- RTMB::sdreport(r$value$obj)$cov.fixed
  expect_identical(unname(vcov(r$value, full = TRUE)), unname(ref))
  cb <- lme4::cbpp
  r <- tier_conds(frm(cbind(incidence, size - incidence) ~ period +
                        (1 | herd), family = binomial(), data = cb))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 0L)
  # a flat sd that costs likelihood when it goes to zero is not a
  # boundary: the at-zero probe says no on a variance the data use
  ns <- asNamespace("frmtmb")
  expect_false(ns$se_at_edge(r$value, "theta_1", -1))
})

test_that("a smooth at its unpenalized limit loses the SE in silence", {
  # mgcv reports a smoothing parameter at infinity without a word, and
  # most GAM fits have one; a smooth is not lme4's singular fit
  set.seed(3)
  d <- data.frame(x = runif(200), z = runif(200))
  d$y <- 2 * d$x + sin(2 * pi * d$z) + rnorm(200, 0, 0.3)
  r <- tier_conds(frm(y ~ s(x) + s(z), data = d))
  lost <- sdr_of(r$value)$se_lost
  # 0.68.1 kept a standard error built on the smooth sd's row of noise
  expect_gt(length(lost), 0L)
  expect_true(all(lost == "boundary"))
  expect_length(r$warnings, 0L)
  expect_length(r$messages, 0L)
})

test_that("an ill-conditioned but identified design keeps every SE", {
  skip_if_not_installed("lme4")
  # the absent case of the eigenvalue gate: a raw polynomial of degree
  # 5 on [1, 2] has a unit-diagonal eigenvalue at 1e-11 of the largest,
  # and the likelihood does curve along it (dev/setier-collin.R)
  set.seed(5)
  d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  fo <- y ~ I(x^1) + I(x^2) + I(x^3) + I(x^4) + I(x^5)
  r <- tier_conds(frm(fo, data = d))
  expect_length(r$warnings, 0L)
  expect_length(sdr_of(r$value)$se_lost, 0L)
  m <- stats::lm(fo, data = d)
  # lm()'s residual variance divides by n - p, the ML one by n
  ref <- sqrt(diag(stats::vcov(m)) * (200 - 6) / 200)
  expect_equal(unname(fixef(r$value)[, "Est.Error"]), unname(ref),
               tolerance = 1e-4)
  # with a random intercept the Hessian is finite differences, and the
  # same holds against lme4
  fo2 <- stats::update(fo, . ~ . + (1 | g))
  r2 <- tier_conds(frm(fo2, data = d))
  expect_length(r2$warnings, 0L)
  expect_length(sdr_of(r2$value)$se_lost, 0L)
  m2 <- allow_warnings(lme4::lmer(fo2, data = d, REML = FALSE),
                       c("unidentifiable", "failed to converge"))
  ref2 <- sqrt(diag(as.matrix(stats::vcov(m2))))
  expect_equal(unname(fixef(r2$value)[, "Est.Error"]), unname(ref2),
               tolerance = 0.01)
})

test_that("a nonlinear ridge with a random effect loses its SEs", {
  # y ~ a + b with a ~ 0 + f: only a_k + b is identified. Without
  # (1 | g) the exact Hessian shows the ridge; with it, tier 1 accepted
  # the finite-difference Hessian, and a prediction of `a` had a
  # standard error of 822,571 (dev/cifixrev-refits.R, R2)
  set.seed(1)
  k <- 10
  d <- data.frame(f = factor(rep(seq_len(k), each = 6)),
                  g = factor(rep(1:6, k)))
  d$y <- rnorm(k)[d$f] + rnorm(6, 0, 0.5)[d$g] + rnorm(k * 6, 0, 0.5)
  fit <- allow_warnings(
    frm(bf(y ~ a + b, a ~ 0 + f, b ~ 1 + (1 | g), nl = TRUE), data = d),
    "are not identified", require = "are not identified")
  lost <- sdr_of(fit)$se_lost
  expect_setequal(names(lost), c(paste0("a_f", seq_len(k)),
                                 "b_(Intercept)"))
  # sigma and the group sd keep theirs
  kept <- !frmtmb:::outer_par_names(fit) %in% names(lost)
  expect_true(all(is.finite(sqrt(diag(sdr_of(fit)$cov.fixed))[kept])))
  w <- character()
  p <- withCallingHandlers(
    frm_linpred(fit, newdata = d[c(1, 7), ], se.fit = TRUE, dpar = "a"),
    frmtmb_se_lost_prediction = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_true(all(is.nan(p$se.fit)))
  expect_length(w, 1L)
})

test_that("separation is named at the default budget, whatever the code", {
  # test-predfix's design, seed 514: z separates yb exactly. nlminb ran
  # to code 0 after 1939 evaluations with the reference BLAS and to code
  # 9 after 1997 with OpenBLAS 0.3.26, where the only warning was "did
  # not report convergence" (dev/reviews/2026-10-06-cifix.md, m5)
  set.seed(514)
  d <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
  d$yb <- as.integer(d$z > 0)
  sep_warn <- function(expr) {
    cls <- character()
    msg <- character()
    v <- withCallingHandlers(expr, warning = function(x) {
      cls <<- c(cls, if (inherits(x, "frmtmb_separation")) "sep" else "other")
      msg <<- c(msg, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
    list(value = v, cls = cls, msg = msg)
  }
  r <- sep_warn(frm(yb ~ z + x, family = bernoulli(), data = d))
  expect_identical(r$cls, "sep")
  expect_match(r$msg, "complete separation", fixed = TRUE)
  expect_match(r$msg, "mu: (Intercept), z", fixed = TRUE)
  # an optimizer stopped early on any platform: the separation is the
  # verdict, not the optimizer's code
  short <- frmtmb_control(optCtrl = list(eval.max = 40, iter.max = 30),
                          restarts = 0)
  r <- sep_warn(frm(yb ~ z + x, family = bernoulli(), data = d,
                    control = short))
  expect_false(r$value$opt$convergence == 0L)
  expect_identical(r$cls, "sep")
  expect_match(r$msg, "where the optimizer stopped (", fixed = TRUE)
})

test_that("quasi-complete separation names its coefficient", {
  # one level of f has no successes: its coefficient runs to -Inf while
  # the others are estimated
  set.seed(8)
  d <- data.frame(f = factor(rep(c("a", "b", "c"), each = 40)),
                  x = rnorm(120))
  d$y <- rbinom(120, 1, plogis(0.3 + 0.5 * d$x))
  d$y[d$f == "c"] <- 0L
  w <- character()
  fit <- withCallingHandlers(frm(y ~ f + x, family = bernoulli(), data = d),
                             frmtmb_separation = function(x) {
                               w <<- c(w, conditionMessage(x))
                               invokeRestart("muffleWarning")
                             })
  expect_length(w, 1L)
  expect_match(w, "quasi-complete separation", fixed = TRUE)
  expect_match(w, "mu: fc", fixed = TRUE)
  # the coefficients the data do hold keep their standard errors
  se <- suppressWarnings(fixef(fit))[, "Est.Error"]
  expect_true(is.finite(se[["x"]]))
})

test_that("a large but finite effect is not separation", {
  # the absent case: a slope of 4 on a covariate spanning about 6 puts
  # fitted probabilities near 1e-5 at both ends, and the maximum exists
  for (s in 1:5) {
    set.seed(s)
    d <- data.frame(x = rnorm(400))
    d$y <- rbinom(400, 1, plogis(-1 + 4 * d$x))
    d$k <- rbinom(400, 5, plogis(1 + 3 * d$x))
    expect_silent(frm(y ~ x, family = bernoulli(), data = d))
    expect_silent(frm(cbind(k, 5 - k) ~ x, family = binomial(), data = d))
  }
})

test_that("ranef(condVar = TRUE) reads the repaired joint covariance", {
  # the gr(g, by = f) fixture of dev/cifixrev-refits.R, seed 11: the
  # slope sd and its correlation lose their standard errors, and
  # sdreport()'s own conditional variances ran from 0 to 99 times what
  # the repaired joint covariance, which predictions read, gives
  set.seed(11)
  dd <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, 10)))
  dd$f <- factor(ifelse(as.integer(dd$g) <= 8, "a", "b"))
  u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
  dd$y <- stats::rnorm(160, 1 + 0.5 * dd$x + u[dd$g, 1] +
                         u[dd$g, 2] * dd$x, 1)
  fit <- suppressMessages(suppressWarnings(
    frm(bf(y ~ x + (1 + x | gr(g, by = f))), family = gaussian(),
        data = dd)))
  skip_if(!length(sdr_of(fit)$se_lost), "no standard error was lost here")
  cv <- as.data.frame(suppressWarnings(ranef(fit, condVar = TRUE)))
  jc <- frmtmb:::get_joint_cov(fit)
  vb <- diag(jc$V)[jc$names == "b"]
  # the long form is ordered by term and coefficient, b by level
  expect_equal(sort(cv$condsd^2), sort(unname(vb)), tolerance = 1e-10)
  # a fit that lost nothing keeps sdreport()'s own conditional variances
  set.seed(1)
  d1 <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  d1$y <- rnorm(100, 1 + 0.5 * d1$x + rnorm(10, 0, 0.8)[d1$g], 1)
  f1 <- frm(y ~ x + (1 | g), data = d1)
  sdr <- sdr_of(f1)
  ref <- sdr$diag.cov.random[names(sdr$par.random) == "b"]
  expect_identical(unname(attr(ranef(f1, condVar = TRUE)$g, "condSD")[, 1]),
                   unname(sqrt(pmax(ref, 0))))
})

# ---- punch round 1 (dev/reviews/2026-10-07-setier.md) -----------------

test_that("a straight step off a curved ridge is not taken as curvature", {
  # B1: c0 + exp(a)^k with a ~ 1 + (1 | g) is invariant to
  # (a0, k, b) -> (s a0, k / s, s b). A straight probe along the ridge's
  # tangent lost 5.6e-3 (above grad_tol) at a ratio of 15.5 per doubling,
  # the fourth power of a curved ridge, and the first probe kept
  # a_(Intercept) at an SE of 9.655 and theta_1 at 0.3677
  set.seed(77)
  dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
  dn$x <- rnorm(nrow(dn))
  u <- rnorm(8, 0, 0.6)[dn$g]
  dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
  f <- suppressMessages(suppressWarnings(
    frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1,
           nl = TRUE), data = dn)))
  lost <- sdr_of(f)$se_lost
  expect_true(all(c("a_(Intercept)", "k_(Intercept)", "theta_1") %in%
                    names(lost)))
  se <- sqrt(diag(vcov(f, full = TRUE)))
  nm <- frmtmb:::outer_par_names(f)
  expect_true(is.nan(se[[match("a_(Intercept)", nm)]]))
  expect_true(is.nan(se[[match("theta_1", nm)]]))
})

test_that("a boundary sd on a larger design is caught by its likelihood", {
  # m1: on 100 groups of 10 the row of an lme4-singular sd rose above the
  # absolute empty-row floor (4.99e-6), and it kept an SE of 1515 with no
  # report (dev/setier-rev-smallsd.R, s1 seed 10)
  skip_if_not_installed("lme4")
  set.seed(10)
  d <- data.frame(g = factor(rep(1:100, each = 10)), x = rnorm(1000))
  d$y <- 1 + 0.5 * d$x + rnorm(100, 0, 0.15)[d$g] + rnorm(1000)
  m <- suppressMessages(lme4::lmer(y ~ x + (1 | g), data = d,
                                   REML = FALSE))
  expect_true(lme4::isSingular(m))
  r <- tier_conds(frm(y ~ x + (1 | g), data = d))
  expect_identical(unname(sdr_of(r$value)$se_lost), "boundary")
  expect_length(r$warnings, 0L)
})

test_that("a code-7 stop at a boundary gets the boundary message", {
  # item 8: y ~ x + (1 + x | g) with no slope variance; nlminb stopped
  # with "singular convergence (7)" at lme4's log-likelihood, every one
  # lme4-singular (dev/setier-rev-rs20.R)
  skip_if_not_installed("lme4")
  set.seed(1)
  d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] + rnorm(120)
  r <- tier_conds(frm(y ~ x + (1 + x | g), data = d))
  expect_false(r$value$opt$convergence == 0L)
  expect_length(r$warnings, 0L)
  bm <- Filter(function(x) inherits(x, "frmtmb_boundary_fit"), r$messages)
  expect_length(bm, 1L)
  # the guard's absent case: with a grad_tol below the gradient left
  # there (1.5e-4) the fit is judged by its gradient and warns
  r2 <- tier_conds(frm(y ~ x + (1 + x | g), data = d,
                       control = frmtmb_control(grad_tol = 1e-12)))
  expect_true(any(grepl("did not report convergence|Large maximum",
                        r2$warnings)))
})

test_that("check_se = 'stop' stops at frm() when the check would wait", {
  # m10: with 22 outer parameters the check waits, and the first
  # summary() stopped instead of frm()
  set.seed(3)
  d <- data.frame(f = factor(rep(1:20, 12)), g = factor(rep(1:12, each = 20)))
  d$y <- rnorm(20)[d$f] + rnorm(240)
  expect_error(frm(y ~ f + (1 | g), data = d,
                   control = frmtmb_control(check_se = "stop")),
               "Boundary (singular) fit", fixed = TRUE)
})

test_that("a separated fit's named coefficients show no SE on any code", {
  # m3: at nlminb's limit the warning said the coefficients have no
  # standard error while vcov() kept sdreport()'s 6.2e132; m5: a probit
  # fit stopped by underflow at a slope of 19.1 was not named
  set.seed(514)
  d <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
  d$yb <- as.integer(d$z > 0)
  short <- frmtmb_control(optCtrl = list(eval.max = 40, iter.max = 30),
                          restarts = 0)
  f <- suppressWarnings(frm(yb ~ z + x, family = bernoulli(), data = d,
                            control = short))
  expect_false(f$opt$convergence == 0L)
  se <- sqrt(diag(vcov(f)))
  expect_true(is.nan(se[["z"]]))
  expect_true(all(sdr_of(f)$se_lost == "separation"))
  set.seed(4)
  dp <- data.frame(x = rnorm(60), z = rnorm(60))
  dp$y <- as.integer(dp$x > 0)
  w <- character()
  withCallingHandlers(frm(y ~ x + z, family = bernoulli("probit"),
                          data = dp),
                      warning = function(x) {
                        w <<- c(w, conditionMessage(x))
                        invokeRestart("muffleWarning")
                      })
  expect_length(w, 1L)
  expect_match(w, "The data separate the outcomes", fixed = TRUE)
})

# ---- punch round 2 (re-check, dev/reviews/2026-10-07-setier.md) -------

test_that("a weak but identified direction keeps its SE off stationarity", {
  # RB1: a poisson pair at a correlation of about 1 - 5e-11. The fit sits
  # 0.07 log-likelihood short along the weak direction, so one probe step
  # gained; a rule on each side's own loss took both SEs, which glm()
  # reports (dev/setier-rev2-window.R, seed 1)
  set.seed(1)
  n <- 400
  d <- data.frame(x1 = rnorm(n))
  d$x2 <- d$x1 + rnorm(n) * 1e-5
  d$y <- rpois(n, exp(-1 + 0.3 * d$x1))
  r <- tier_conds(frm(y ~ x1 + x2, family = poisson(), data = d))
  expect_length(sdr_of(r$value)$se_lost, 0L)
  expect_length(r$warnings, 0L)
  ref <- sqrt(diag(stats::vcov(stats::glm(y ~ x1 + x2, family = poisson(),
                                          data = d))))
  se <- fixef(r$value)[, "Est.Error"]
  # base frmtmb and glm() differ by up to 2 percent here (both invert a
  # Hessian conditioned to 1e10)
  expect_equal(unname(se), unname(ref), tolerance = 0.05)
})

test_that("only a proved separation relabels a flat coefficient", {
  # m-c: a flat coefficient above 10 in a bernoulli polynomial on data
  # that are not separated was called "separation"
  ns <- asNamespace("frmtmb")
  set.seed(1)
  d <- data.frame(x = runif(400, 1, 2))
  d$y <- rbinom(400, 1, plogis(-1 + sin(3 * d$x)))
  f <- frm(y ~ x, family = bernoulli(), data = d)
  nm <- ns$outer_par_names(f)
  lost <- stats::setNames("flat", nm[2])
  f$opt$par[2] <- 50
  expect_identical(unname(ns$se_relabel(f, lost)), "flat")
  f$cache$se_explained_sep <- nm
  expect_identical(unname(ns$se_relabel(f, lost)), "separation")
})

test_that("an sd the fit left near zero short of the maximum is no boundary", {
  # RB2: with the response in units of 1e-3, seed 36, nlminb stopped with
  # "false convergence (8)" at an sd of 1.7e-12 times sigma, 0.825 below
  # lme4's log-likelihood at 0.33; the convergence warning was replaced
  # by "is at zero" (dev/setier-rev2-trap.R)
  set.seed(36)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- (1 + 0.5 * d$x + rnorm(100)) * 1e-3
  r <- tier_conds(frm(y ~ x + (1 | g), data = d))
  bm <- Filter(function(x) inherits(x, "frmtmb_boundary_fit"), r$messages)
  expect_length(bm, 0L)
  expect_true(any(grepl("did not report convergence|stopped short",
                        r$warnings)))
  # seed 19, code 0, 0.666 short: the SE warning says why
  set.seed(19)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- (1 + 0.5 * d$x + rnorm(100)) * 1e-3
  r <- tier_conds(frm(y ~ x + (1 | g), data = d))
  expect_length(Filter(function(x) inherits(x, "frmtmb_boundary_fit"),
                       r$messages), 0L)
  expect_true(any(grepl("stopped short of the maximum", r$warnings,
                        fixed = TRUE)))
})

test_that("under quadrature no sd is called a boundary", {
  # final check RQ1: the upward check cannot run on a quadrature
  # objective, and seed 36 at x1e-3 refitted with quadrature = TRUE,
  # 0.825 below lme4's log-likelihood, was told "is at zero"
  set.seed(36)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- (1 + 0.5 * d$x + rnorm(100)) * 1e-3
  r <- tier_conds(frm(y ~ x + (1 | g), data = d, quadrature = TRUE))
  bm <- Filter(function(x) inherits(x, "frmtmb_boundary_fit"), r$messages)
  expect_length(bm, 0L)
  expect_false(any(sdr_of(r$value)$se_lost %in% c("boundary", "short")))
})
