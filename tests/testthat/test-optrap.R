# Fits that stopped short of the maximum likelihood, or ended on luck,
# without saying so (lane optrap, dev/optrap-findings.md).

# ---- the units of the response --------------------------------------

# y ~ x + (1 | g) with no group effect in the truth, 20 groups of 5:
# the design on which a rescaled response stopped short (lane setier's
# review, RB2). Seeds 19, 36 and 37 have an interior lme4 maximum.
optrap_ri20 <- function(seed) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(100)
  d
}

# sleepstudy's shape: 18 subjects, 10 days, a random slope
optrap_slope <- function(seed) {
  set.seed(seed)
  d <- data.frame(s = factor(rep(1:18, each = 10)), t = rep(0:9, 18))
  b0 <- rnorm(18, 0, 25)
  b1 <- rnorm(18, 0, 6)
  d$y <- (250 + b0[d$s] + (10 + b1[d$s]) * d$t + rnorm(180, 0, 25)) / 25
  d
}

# The fit of y * k must be the fit of y with every coefficient of mu, the
# residual sd and the group sds multiplied by k: the density of y * k is
# the density of y over k^n. Compared on the scale-1 fit's own
# yardsticks: each estimate within 1e-2 of its standard error, and the
# log-likelihood within what moving the group log sd by 1e-2 of its
# standard error costs at the scale-1 optimum.
optrap_invariance <- function(form, d, scales = 10^c(-6, -3, 3, 6)) {
  n <- nrow(d)
  pick <- function(f) {
    c(f$estimates$beta, f$estimates$betad[1], f$estimates$theta[1])
  }
  shift <- function(k) {
    c(rep(k, 2), log(k), log(k))
  }
  f1 <- frm(form, data = d)
  # on the scale the optimizer iterates on (log sd, log sigma)
  H <- stats::optimHess(f1$opt$par, f1$obj$fn, f1$obj$gr)
  se1 <- sqrt(diag(solve(H)))
  nm <- frmtmb:::outer_par_names(f1)
  j <- match("theta_1", nm)
  se <- c(se1[seq_len(2)], se1[match("sigma_(Intercept)", nm)], se1[j])
  expect_true(all(is.finite(se)))
  p <- f1$opt$par
  p[j] <- p[j] + 0.01 * se[4]
  tol_ll <- f1$obj$fn(p) - f1$opt$objective
  expect_gt(tol_ll, 0)
  ll1 <- as.numeric(logLik(f1))
  est1 <- pick(f1)
  for (k in scales) {
    dk <- d
    dk$y <- d$y * k
    fk <- frm(form, data = dk)
    expect_identical(fk$opt$convergence, 0L)
    expect_lt(abs(as.numeric(logLik(fk)) + n * log(k) - ll1), tol_ll)
    ek <- pick(fk)
    sk <- shift(k)
    back <- c(ek[1:2] / sk[1:2], ek[3:4] - sk[3:4])
    expect_lt(max(abs(back - est1) / se), 1e-2)
    if (requireNamespace("lme4", quietly = TRUE)) {
      l4 <- suppressMessages(lme4::lmer(form, data = dk, REML = FALSE))
      expect_gt(as.numeric(logLik(fk)) - as.numeric(logLik(l4)), -tol_ll)
    }
  }
}

test_that("a random intercept fit is the same at every scale of y", {
  skip_on_cran()
  # rellib-r7: at y x 1e3 seeds 19, 36 and 37 stopped 0.67, 0.82 and
  # 1.68 log-likelihood units short of lme4 at code 0
  for (s in c(19, 36, 37)) {
    optrap_invariance(y ~ x + (1 | g), optrap_ri20(s))
  }
})

test_that("a random slope fit is the same at every scale of y", {
  skip_on_cran()
  # rellib-r7: 17 of 40 such fits at y x 1e3 and 40 of 40 at 1e6 short
  optrap_invariance(y ~ t + (1 + t | s), optrap_slope(2))
})

test_that("only a predictor in the response's units gets a unit", {
  d <- optrap_ri20(19)
  fr <- function(form, family, data) {
    frm(form, family = family, data = data, dry_run = "frame")
  }
  units <- function(f) frmtmb:::predictor_units(f)
  # a response near unit spread keeps the fit it always had
  expect_length(units(fr(y ~ x + (1 | g), gaussian(), d)), 0L)
  d3 <- d
  d3$y <- d$y * 1e3
  u <- units(fr(y ~ x + (1 | g), gaussian(), d3))
  expect_length(u, 1L)
  # a power of two, within a factor sqrt(2) of the spread the fixed
  # design leaves (predictor_resid_sd())
  expect_identical(log2(u[[1]]), round(log2(u[[1]])))
  r <- stats::residuals(stats::lm(y ~ x, data = d3))
  expect_lt(abs(log2(u[[1]] / stats::mad(r))), 0.5 + 1e-12)
  # a log-link count and a lognormal's log scale are in no such units
  dp <- d
  dp$y <- rpois(100, 3) * 1000
  expect_length(units(fr(y ~ x + (1 | g), poisson(), dp)), 0L)
  dl <- d
  dl$y <- exp(d$y) * 1e3
  expect_length(units(fr(y ~ x + (1 | g), lognormal(), dl)), 0L)
  # sigma's own predictor is on its log link
  uf <- units(fr(bf(y ~ x + (1 | g), sigma ~ x), gaussian(), d3))
  expect_identical(names(uf), names(u))
  # a sigma held by se() is no estimated scale to take a unit from
  d3$s <- 1e3
  expect_length(units(fr(y | se(s) ~ x + (1 | g), gaussian(), d3)), 0L)
  expect_length(units(fr(y | se(s) ~ x + (1 | g), student(), d3)), 0L)
  # an exact fit leaves only rounding, which is no scale: a constant
  # response, and one the design reproduces
  dc <- d
  dc$y <- 2
  expect_length(units(fr(y ~ x + (1 | g), gaussian(), dc)), 0L)
  dc$y <- 3 + 1e3 * d$x
  expect_length(units(fr(y ~ x + (1 | g), gaussian(), dc)), 0L)
})

test_that("a covariate that explains most of y leaves the fit silent", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  # y = 1 + 1e4 x + b_g + e: sd(y) is near 1e4 and the residual sd near
  # 1. A unit from sd(y), lane optrap's first rule, left these fits at
  # lme4's maximum but at "false convergence (8)" with a warning, 34 of
  # 40 seeds, these three among them (dev/optrap-r2.R with
  # OPTRAP_ROUND0=1; the review of lane optrap, B1)
  for (s in 1:3) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
    d$y <- 1 + 1e4 * d$x + rnorm(20, 0, 0.5)[d$g] + rnorm(100)
    w <- character()
    f <- withCallingHandlers(frm(y ~ x + (1 | g), data = d),
                             warning = function(x) {
                               w <<- c(w, conditionMessage(x))
                               invokeRestart("muffleWarning")
                             })
    expect_length(w, 0L)
    expect_identical(f$opt$convergence, 0L)
    # at lme4's maximum, within what moving the residual log sd by 1e-2
    # of its standard error costs here
    H <- stats::optimHess(f$opt$par, f$obj$fn, f$obj$gr)
    j <- match("betad", names(f$opt$par))
    p <- f$opt$par
    p[j] <- p[j] + 0.01 * sqrt(solve(H)[j, j])
    tol_ll <- f$obj$fn(p) - f$opt$objective
    l4 <- suppressMessages(lme4::lmer(y ~ x + (1 | g), data = d,
                                      REML = FALSE))
    expect_gt(as.numeric(logLik(f)) - as.numeric(logLik(l4)), -tol_ll)
  }
})

test_that("a few outliers do not set the unit", {
  # t(1) residuals, seed 18 of the review's heavy-tail case: the sd of
  # the residuals gave a unit of 64 against a scale near 1, and the
  # student fit stopped 29.8 short (the review of lane optrap, r1); the
  # MAD gives the scale of the bulk
  set.seed(18)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.5)[d$g] + rt(100, 1)
  fr <- frm(y ~ x + (1 | g), family = student(), data = d,
            dry_run = "frame")
  u <- frmtmb:::predictor_units(fr)
  r <- stats::residuals(stats::lm(y ~ x, data = d))
  expect_lt(abs(log2((if (length(u)) u[[1]] else 1) / stats::mad(r))),
            0.5 + 1e-12)
  expect_gt(stats::sd(r) / stats::mad(r), 8)
})

test_that("a design above the least-squares size cap keeps the residual unit", {
  # Above the cap the least squares runs on evenly spaced rows. Punch
  # round 1 skipped it there and took the raw spread of y: 101000 rows
  # by 50 columns got a unit of 8192 where the residual spread is 1, and
  # ran to the evaluation limit 8.5e3 short (the review of lane optrap,
  # dev/optrap-rev2-b1.R wide). The cap is lowered here so that a small
  # design crosses it: 200 cells against 150, so 75 rows of 100 set the
  # spread, and the unit is the full design's within one power of two
  # (the raw spread would give 8192)
  unit1 <- function(fr) {
    u <- frmtmb:::predictor_units(fr)
    if (length(u)) u[[1]] else 1
  }
  for (s in 1:3) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
    d$y <- 1 + 1e4 * d$x + rnorm(20, 0, 0.5)[d$g] + rnorm(100)
    fr <- frm(y ~ x + (1 | g), data = d, dry_run = "frame")
    full <- unit1(fr)
    old <- options(frmtmb.unit_ls_cells = 150)
    capped <- unit1(fr)
    options(old)
    expect_lte(abs(log2(capped / full)), 1)
  }
})

test_that("an ill-conditioned fixed-effect fit reaches lm()'s maximum", {
  # A raw polynomial of degree 6 on x in [1, 2]: nlminb's relative
  # function test passed 0.5 to 1 standard errors from lm() in every
  # coefficient, at code 0 and in silence, on 19 of 20 seeds with the
  # residual unit and on 19 of 20 on 0.69.0 with y times 16
  # (dev/optrap-p1-polystudy.R); newton_polish() takes the rest
  for (s in 1:2) for (k in c(1, 16)) {
    set.seed(6000 + s)
    d <- data.frame(x = runif(200, 1, 2))
    d$y <- (sin(3 * d$x) + rnorm(200, 0, 0.2)) * k
    for (j in 1:6) d[[paste0("x", j)]] <- d$x^j
    fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6
    w <- character()
    f <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
    expect_length(w, 0L)
    expect_identical(f$opt$convergence, 0L)
    l <- stats::lm(fo, data = d)
    se <- sqrt(diag(stats::vcov(l)))
    expect_lt(max(abs(unname(fixef_by_dpar(f)$mu) - unname(coef(l))) /
                    se), 0.01)
  }
})

test_that("a Newton step needs a positive definite Hessian", {
  # f = x^2 - y^2 + 0.1 y has no minimum; from (1, 0.1) the Newton step
  # goes to the saddle (0, 0.05) and f falls on the way (the review of
  # lane optrap, r3)
  obj <- RTMB::MakeADFun(function(p) p$x^2 - p$y^2 + 0.1 * p$y,
                         list(x = 1, y = 0.1), silent = TRUE)
  opt <- list(par = c(x = 1, y = 0.1), objective = obj$fn(c(1, 0.1)),
              convergence = 0L)
  out <- frmtmb:::newton_polish(obj, opt, list(lower = -Inf, upper = Inf))
  expect_identical(out, opt)
})

test_that("a caller's optimizer tolerance leaves the fit where it stopped", {
  # optim's own names count as nlminb's do (the review of lane optrap,
  # r4); the design is the one newton_polish() steps on by default
  set.seed(6001)
  d <- data.frame(x = runif(200, 1, 2))
  d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
  for (j in 1:6) d[[paste0("x", j)]] <- d$x^j
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6
  for (ctl in list(list(reltol = 1e-8), list(factr = 1e7))) {
    f <- suppressWarnings(frm(fo, data = d, control = frmtmb_control(
      optimizer = "optim", optCtrl = ctl)))
    expect_null(f$opt$polish)
  }
  f0 <- suppressWarnings(frm(fo, data = d, control = frmtmb_control(
    optimizer = "optim")))
  expect_false(is.null(f0$opt$polish))
})

test_that("a Newton step that gains far less than it predicts is refused", {
  # negative binomial data at the Poisson limit, where the shape's
  # maximum is at infinity; seed 2 of the 0.70.0 consolidation's
  # dev/rel070-polishgate.R, with its draws in its order. A step that
  # gained 3e-8 of a prediction above 1e-6 took the log shape from 20 to
  # 67, where the gradient is 565, and both coefficients lost their SEs
  set.seed(2)
  invisible(runif(200)); invisible(rnorm(20)); invisible(rnorm(200))
  set.seed(2)
  n <- 150
  x <- rnorm(n)
  invisible(rnorm(n))
  mu <- stats::plogis(-2 + 0.5 * x)
  invisible(rbeta(n, mu * 200, (1 - mu) * 200)); invisible(sample(n, 3))
  invisible(rbeta(n, 0.98 * 20, 0.02 * 20)); invisible(sample(n, 2))
  invisible(rt(n, df = 30))
  dn <- data.frame(x = x, c = rpois(n, exp(1 + 0.3 * x)))
  fit <- suppressWarnings(suppressMessages(
    frm(c ~ x, family = negbinomial(), data = dn)))
  expect_null(fit$opt$polish)
  se <- suppressWarnings(fixef(fit)[, "Est.Error"])
  expect_true(all(is.finite(se[1:2])))
})

test_that("Newton steps on an ill-conditioned design reach the maximum", {
  # lane seflags' B5 design without its group term, seeds 4 and 5:
  # degree-8 raw powers plus z and w. In natural units the steps gained
  # what they predicted and ended where the Hessian has no Cholesky
  # factor (a gradient of 0.0335 and 1.63 on sigma), and were undone, so
  # the fits stayed 0.674 and 1.73 below the maximum (dev/optrap-c3-b5.R,
  # filed at 0.70.0). Since lane illse they are taken in the design's
  # coordinates, where the Hessian factors (reparam_basis()), and reach
  # the maximum of the same likelihood with poly(x, 8)
  for (s in c(4, 5)) {
    set.seed(s)
    d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)),
                    z = stats::rnorm(200),
                    w = factor(sample(c("a", "b"), 200, TRUE)))
    d$y <- sin(3 * d$x) + 0.3 * d$z + 0.5 * (d$w == "b") +
      stats::rnorm(20, 0, 0.3)[d$g] + stats::rnorm(200, 0, 0.2)
    for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
    w <- character()
    f <- withCallingHandlers(
      frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + z + w, data = d),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    expect_false(is.null(f$opt$polish))
    expect_false(any(grepl("Large maximum absolute gradient", w,
                           fixed = TRUE)))
    ref <- frm(y ~ poly(x, 8) + z + w, data = d)
    # within a tenth of what one SE of log sigma moves the likelihood
    expect_lt(abs(as.numeric(logLik(ref) - logLik(f))), 0.05)
  }
})

test_that("a Newton step the quadratic model overrates is still taken", {
  # degree-7 raw powers, seed 17 of dev/optrap-c3-repoly.R, no group
  # term: the first step gains 0.2 of its prediction (log sigma is not a
  # quadratic direction) and the later ones 1.00 to 1.06; a gate of a
  # quarter refused the first and left the fit 1.16 short in silence
  # (dev/optrap-c3-trace17.R)
  set.seed(7017)
  d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7
  f <- suppressWarnings(frm(fo, data = d))
  expect_false(is.null(f$opt$polish))
  l <- stats::lm(fo, data = d)
  se <- sqrt(diag(stats::vcov(l)))
  expect_lt(max(abs(unname(fixef_by_dpar(f)$mu) - unname(coef(l))) / se),
            0.01)
})

test_that("a near-collinear fit with a group term reaches the maximum", {
  skip_if_not_installed("lme4")
  # raw powers of x on [1, 2] to degree 6 with (1 | g): nlminb stopped
  # short at code 0 with no warning on 13 of 20 seeds, these up to 1.1
  # (dev/optrap-c3-repoly.R); lme4 reaches the maximum on all 20
  for (s in c(1, 13, 14)) {
    set.seed(6000 + s)
    d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
    d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
    for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
    fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + (1 | g)
    w <- character()
    f <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
    expect_length(w, 0L)
    # within what moving the residual log sd by 1e-2 of its standard
    # error costs here
    H <- stats::optimHess(f$opt$par, f$obj$fn, f$obj$gr,
                          control = list(parscale = f$par_units %||%
                                           rep(1, length(f$opt$par))))
    j <- match("betad", names(f$opt$par))
    p <- f$opt$par
    p[j] <- p[j] + 0.01 * sqrt(solve(H)[j, j])
    tol_ll <- f$obj$fn(p) - f$opt$objective
    l4 <- suppressMessages(lme4::lmer(fo, data = d, REML = FALSE))
    expect_gt(as.numeric(logLik(f)) - as.numeric(logLik(l4)), -tol_ll)
  }
})

# ---- a restart is a second chance, not a requirement -----------------

# A taped objective with a minimum at a = 3, b = -1, and a gradient that
# raises "NA/NaN gradient evaluation" once `fail_after` evaluations have
# been made, as RTMB does where a step leaves the likelihood's domain.
optrap_obj <- function() {
  obj <- RTMB::MakeADFun(function(p) {
    (p$a - 3)^2 + 10 * (p$b + 1 - 0.1 * p$a^2)^2
  }, list(a = 0, b = 0), silent = TRUE)
  n <- 0L
  fail_after <- Inf
  bump <- 0
  w <- obj
  w$fn <- function(x) {
    n <<- n + 1L
    obj$fn(x) + if (n > fail_after) bump else 0
  }
  w$gr <- function(x) {
    n <<- n + 1L
    if (n > fail_after && bump == 0) stop("NA/NaN gradient evaluation")
    obj$gr(x)
  }
  w$set <- function(after, b = 0) {
    fail_after <<- after
    bump <<- b
  }
  w$count <- function() n
  w
}

test_that("a restart that raises keeps the run before it", {
  # rellib-r7: the error left optimize_obj() and took the fit with it,
  # although the first run's optimum was finite (fuzz seed 20379118
  # under OpenBLAS's LAPACK)
  ctl <- frmtmb_control(optCtrl = list(iter.max = 3, eval.max = 100))
  first <- optrap_obj()
  op1 <- frmtmb:::optimize_obj(first, utils::modifyList(ctl,
                                                          list(restarts = 0)))
  # the case needs a first run that stops above grad_tol
  expect_gt(max(abs(first$gr(op1$par))), ctl$grad_tol)
  w <- optrap_obj()
  w$set(first$count())
  op <- frmtmb:::optimize_obj(w, ctl)
  expect_identical(unname(op$par), unname(op1$par))
  expect_identical(op$objective, op1$objective)
  lf <- w$env$lfixed()
  expect_identical(unname(w$env$last.par[lf]), unname(op$par))
  expect_identical(unname(w$env$last.par.best[lf]), unname(op$par))
})

test_that("a restart that ends worse leaves the objective at the kept run", {
  # rellib-r7: the restart's own last point stayed the objective's, so
  # parList() would read the random effects' modes there
  ctl <- frmtmb_control(optCtrl = list(iter.max = 3, eval.max = 100))
  first <- optrap_obj()
  op1 <- frmtmb:::optimize_obj(first, utils::modifyList(ctl,
                                                          list(restarts = 0)))
  w <- optrap_obj()
  # every evaluation of the restart costs 100 more, so it ends worse
  w$set(first$count(), b = 100)
  op <- frmtmb:::optimize_obj(w, ctl)
  expect_identical(unname(op$par), unname(op1$par))
  lf <- w$env$lfixed()
  expect_identical(unname(w$env$last.par[lf]), unname(op$par))
  expect_identical(unname(w$env$last.par.best[lf]), unname(op$par))
  # the inverse: a restart that improves is kept, and is where it stopped
  w2 <- optrap_obj()
  op2 <- frmtmb:::optimize_obj(w2, ctl)
  expect_lt(op2$objective, op1$objective)
  expect_identical(unname(w2$env$last.par[lf]), unname(op2$par))
})

# ---- mo(): the other sign's cone ----------------------------------------

# brms_monotonic's data code and the exact maximum of
# `ls ~ mo(income) * age` over every sign and support, as in
# test-optima.R
optrap_mo_data <- function(seed) {
  set.seed(seed)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}

optrap_mo_exact <- function(d) {
  code <- as.integer(d$income) - 1L
  S <- sapply(1:3, function(k) as.numeric(code >= k))
  n <- nrow(d)
  best <- -Inf
  for (s1 in c(-1, 1)) for (s2 in c(-1, 1)) for (m1 in 0:7) for (m2 in 0:7) {
    a1 <- which(bitwAnd(m1, c(1L, 2L, 4L)) > 0)
    a2 <- which(bitwAnd(m2, c(1L, 2L, 4L)) > 0)
    X <- cbind(1, d$age, S[, a1, drop = FALSE],
               (S * d$age)[, a2, drop = FALSE])
    f <- stats::lm.fit(X, d$ls)
    cf <- f$coefficients[-(1:2)]
    if (anyNA(cf)) next
    sg <- c(rep(s1, length(a1)), rep(s2, length(a2)))
    if (any(sg * cf < 0)) next
    best <- max(best,
                -(n / 2) * (log(2 * pi * sum(f$residuals^2) / n) + 1))
  }
  best
}

test_that("the other cone's search reaches its maximum", {
  # rellib-r7: seed 194 stopped 0.0489 and seed 354 0.449 below the
  # maximum at code 0 and no warning; seeds 80 (a coefficient near 0)
  # and 184 stopped short on builds of this lane before the search
  # covered them (dev/optrap-findings.md)
  for (s in c(80, 184, 194, 354)) {
    d <- optrap_mo_data(s)
    fit <- frm(ls ~ mo(income) * age, data = d)
    ex <- optrap_mo_exact(d)
    expect_identical(fit$opt$convergence, 0L)
    expect_lt((ex - as.numeric(logLik(fit))) / abs(ex), 1e-8,
              label = paste("relative gap below the maximum, seed", s))
  }
})
