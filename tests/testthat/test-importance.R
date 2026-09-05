# Importance-sampling correction of the Laplace approximation
# (frm(importance =), R/importance.R). Reference validation against
# GLMMadaptive and lme4 runs locally and in CI, not on CRAN.
skip_on_cran()

# The probe design the correction was built for: 60 groups of 8
# Bernoulli rows with a correlated random intercept and slope, which is
# where the Laplace approximation is measurably biased and where
# quadrature = TRUE cannot go (its transform is one-dimensional).
imp_probe_data <- function() {
  set.seed(11)
  ng <- 60L
  per <- 8L
  g <- rep(seq_len(ng), each = per)
  x <- rnorm(ng * per)
  sig <- matrix(c(1, 0.35, 0.35, 0.49), 2)
  b <- MASS::mvrnorm(ng, c(0, 0), sig)
  data.frame(y = rbinom(ng * per, 1,
                        plogis(-0.5 + 0.8 * x + b[g, 1] + b[g, 2] * x)),
             x = x, g = factor(g))
}

# A gaussian design, where the Laplace approximation is EXACT.
imp_gauss_data <- function() {
  set.seed(5)
  ng <- 50L
  per <- 6L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u <- rnorm(ng, 0, 0.8)
  data.frame(y = rnorm(ng * per, 1 + 0.5 * x + u[g], 1), x = x, g = g)
}

# The scalar-intercept design the quadrature tests already use.
imp_scalar_data <- function() {
  set.seed(601)
  ng <- 100L
  per <- 4L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u <- rnorm(ng, 0, 1.2)
  data.frame(y = rbinom(ng * per, 1, plogis(-0.5 + 0.7 * x + u[g])),
             x = x, g = g)
}

# The Laplace objective, its optimum, and the pieces the correction is
# built from, without going through frm(). The unit tests below poke at
# the proposal directly, which is the only way to hold the anchor still
# while the parameters move.
imp_parts <- function(bform, data) {
  fr <- frm(bform, data = data, dry_run = "frame")
  nll <- frmtmb:::build_objective(fr)
  tpl <- frmtmb:::make_start(fr, NULL, NULL)
  lap <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr[["map"]],
                         silent = TRUE)
  opt <- nlminb(lap$par, lap$fn, lap$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  bk <- fr[["re_blocks"]][[1L]]
  list(frame = fr, nll = nll, tpl = tpl, lap = lap, opt = opt, bk = bk,
       gmap = frmtmb:::imp_group_map(fr, bk))
}

# ---------------------------------------------------------------------
# (a) THE GAUSSIAN IDENTITY, the free pin of the plumbing.
#
# For a gaussian response the Laplace approximation is exact, so every
# importance weight is equal and the correction has nothing to correct.
# The identity holds AT the proposal's anchor, which is where the fit
# reports from, and it holds for ANY draw count because it is an
# algebraic cancellation and not a limit.
# ---------------------------------------------------------------------
test_that("gaussian: the corrected objective at its anchor IS Laplace's", {
  p <- imp_parts(bf(y ~ x + (1 | g)) + gaussian(), imp_gauss_data())
  for (nd in c(2L, 10L, 1000L)) {
    plan <- frmtmb:::imp_plan(p$lap, p$frame, p$bk, p$opt$par, nd, 1L)
    io <- frmtmb:::build_importance_objective(p$frame, p$bk, p$gmap, plan)
    pl <- frmtmb:::imp_par_list(p$tpl, p$opt$par)
    expect_equal(io$fn(pl), p$opt$objective, tolerance = 1e-10)
    # every weight equal, so the effective sample size is the whole
    # sample and the Monte Carlo error is exactly zero
    ess <- frmtmb:::imp_ess(as.matrix(io$amat(pl)), plan[["n_draw"]])
    expect_equal(min(ess$ess), 1, tolerance = 1e-8)
    expect_equal(ess$mcse, 0, tolerance = 1e-8)
  }
})

test_that("gaussian: the fitted correction reproduces the Laplace fit", {
  dd <- imp_gauss_data()
  fl <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  fi <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
            importance = 1000L)
  # Laplace is exact here, so the correction must not move the answer.
  # What is left is the optimizer's own drift on a Monte Carlo surface.
  expect_lt(abs(as.numeric(logLik(fi)) - as.numeric(logLik(fl))), 1e-4)
  expect_lt(max(abs(fi$opt$par - fl$opt$par)), 1e-3)
  expect_equal(fi$importance$ess_min, 1, tolerance = 1e-6)
})

# ---------------------------------------------------------------------
# (b) THE PROBE DESIGN. The Laplace fit is biased; the correction must
# land on the adaptive Gauss-Hermite answer within the Monte Carlo
# standard error the fit itself reports.
# ---------------------------------------------------------------------
test_that("the probe design lands on GLMMadaptive within its MCSE", {
  skip_if_not_installed("GLMMadaptive")
  skip_if_not_installed("MASS")
  dd <- imp_probe_data()
  ga <- GLMMadaptive::mixed_model(y ~ x, ~ x | g, dd, binomial(),
                                  nAGQ = 25)
  ref_nll <- -as.numeric(logLik(ga))
  ref_fix <- unname(GLMMadaptive::fixef(ga))

  fl <- frm(bf(y ~ x + (x | g)) + binomial(), data = dd)
  fi <- frm(bf(y ~ x + (x | g)) + binomial(), data = dd,
            importance = 2000L, se = TRUE)
  mcse <- fi$importance$mcse
  expect_gt(mcse, 0)

  # the corrected log-likelihood, against its own reported MCSE
  expect_lt(abs(fi$opt$objective - ref_nll), 3 * mcse)
  # the intercept, against its own standard error
  se_int <- sqrt(diag(vcov(fi)))[[1L]]
  expect_lt(abs(fixef(fi)$mu[[1L]] - ref_fix[[1L]]), 0.5 * se_int)
  # the Laplace fit it corrects is genuinely biased: it misses the
  # reference log-likelihood by more than a unit, far outside the MCSE
  expect_gt(abs(fl$opt$objective - ref_nll), 1)
  expect_gt(abs(fl$opt$objective - ref_nll), 20 * mcse)
  # and the correction closes almost all of that gap
  expect_lt(abs(fi$opt$objective - ref_nll),
            0.05 * abs(fl$opt$objective - ref_nll))
})

# ---------------------------------------------------------------------
# (c) SCALAR INTERCEPT: agreement with this package's own quadrature
# and with lme4's adaptive Gauss-Hermite, on the design the suite
# already uses for that comparison.
# ---------------------------------------------------------------------
test_that("scalar intercept agrees with quadrature and glmer(nAGQ = 25)", {
  skip_if_not_installed("lme4")
  dd <- imp_scalar_data()
  fq <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd,
            quadrature = TRUE)
  fi <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd,
            importance = 1000L)
  ref <- lme4::glmer(y ~ x + (1 | g), dd, family = binomial, nAGQ = 25)

  mcse <- fi$importance$mcse
  expect_lt(abs(as.numeric(logLik(fi)) - as.numeric(logLik(ref))),
            3 * mcse)
  expect_lt(abs(as.numeric(logLik(fi)) - as.numeric(logLik(fq))),
            3 * mcse)
  expect_lt(max(abs(fixef(fi)$mu - lme4::fixef(ref))), 0.05)
  sd_i <- sqrt(VarCorr(fi)[[1L]][1, 1])
  sd_r <- as.numeric(attr(lme4::VarCorr(ref)$g, "stddev"))[1L]
  expect_lt(abs(sd_i - sd_r), 0.05)
  # the conditional modes survive: the corrected tape carries none, so
  # they come from the Laplace inner solve at the corrected optimum
  re <- ranef(fi)[[1L]]
  expect_false(anyNA(unlist(re)))
  expect_identical(nrow(re), 100L)
})

# ---------------------------------------------------------------------
# (d) THE BIAS DEMONSTRATION. Binary data in small clusters, where the
# Laplace approximation shrinks the variance component and the
# correction recovers the generating value.
# ---------------------------------------------------------------------
test_that("the correction recovers a variance component Laplace shrinks", {
  set.seed(404)
  ng <- 120L
  per <- 3L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  sd_true <- 1.5
  u <- rnorm(ng, 0, sd_true)
  dd <- data.frame(y = rbinom(ng * per, 1, plogis(-0.3 + 0.6 * x + u[g])),
                   x = x, g = g)
  fl <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd)
  fi <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd,
            importance = 2000L)
  sd_lap <- sqrt(VarCorr(fl)[[1L]][1, 1])
  sd_imp <- sqrt(VarCorr(fi)[[1L]][1, 1])
  # Laplace is biased DOWN here, and the correction moves up toward the
  # generating value
  expect_lt(sd_lap, sd_true)
  expect_gt(sd_imp, sd_lap)
  expect_lt(abs(sd_imp - sd_true), abs(sd_lap - sd_true))
  # the corrected log-likelihood is higher than the Laplace one by far
  # more than the Monte Carlo error of the correction
  expect_lt(fi$opt$objective, fl$opt$objective)
  expect_gt(fl$opt$objective - fi$opt$objective, 5 * fi$importance$mcse)
})

# ---------------------------------------------------------------------
# (e) THE GRADIENT, against numDeriv, at the optimum and away from it.
# ---------------------------------------------------------------------
test_that("the corrected objective's gradient matches numDeriv", {
  skip_if_not_installed("numDeriv")
  p <- imp_parts(bf(y ~ x + (1 | g)) + binomial(), imp_scalar_data())
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$bk, p$opt$par, 200L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$bk, p$gmap, plan)
  otpl <- frmtmb:::imp_template(p$tpl, "b")
  obj <- RTMB::MakeADFun(io$fn, otpl, silent = TRUE)
  for (shift in c(0, 0.3)) {
    at <- p$opt$par + shift
    nd <- numDeriv::grad(obj$fn, at)
    expect_equal(as.numeric(obj$gr(at)), nd, tolerance = 1e-6)
  }
})

# ---------------------------------------------------------------------
# (f) THE ESS DIAGNOSTIC fires on a proposal frozen deliberately far
# from the optimum and stays quiet at the optimum.
# ---------------------------------------------------------------------
test_that("the effective sample size diagnostic separates the two regimes", {
  skip_if_not_installed("MASS")
  p <- imp_parts(bf(y ~ x + (x | g)) + binomial(), imp_probe_data())
  ess_at <- function(anchor, at) {
    plan <- frmtmb:::imp_plan(p$lap, p$frame, p$bk, anchor, 1000L, 1L)
    io <- frmtmb:::build_importance_objective(p$frame, p$bk, p$gmap, plan)
    frmtmb:::imp_ess(as.matrix(io$amat(frmtmb:::imp_par_list(p$tpl, at))),
                     plan[["n_draw"]])
  }
  at_opt <- ess_at(p$opt$par, p$opt$par)
  expect_gt(min(at_opt$ess), frmtmb:::imp_ess_floor)
  # displace the proposal's WIDTHS by half a log standard deviation.
  # The estimate stays unbiased, so what degrades is efficiency, and
  # that is exactly what this diagnostic is for.
  far <- p$opt$par
  far[3:4] <- far[3:4] - 0.5
  at_far <- ess_at(far, p$opt$par)
  expect_lt(min(at_far$ess), frmtmb:::imp_ess_floor)
  expect_gt(at_far$mcse, 5 * at_opt$mcse)
  # and the warning names the groups
  expect_warning(
    frmtmb:::imp_ess_warning(at_far$ess, p$bk, 1000L,
                             frmtmb:::imp_ess_floor),
    "covers"
  )
  expect_silent(
    frmtmb:::imp_ess_warning(at_opt$ess, p$bk, 1000L,
                             frmtmb:::imp_ess_floor)
  )
})

# ---------------------------------------------------------------------
# (g) EVERY REFUSAL, by name.
# ---------------------------------------------------------------------
test_that("importance refuses the arguments it cannot honor", {
  dd <- imp_gauss_data()
  bfm <- bf(y ~ x + (1 | g)) + gaussian()
  expect_error(frm(bfm, dd, importance = -1), "importance")
  expect_error(frm(bfm, dd, importance = 2.5), "importance")
  expect_error(frm(bfm, dd, importance = 100L, quadrature = TRUE),
               "two different")
  expect_error(frm(bfm, dd, importance = 100L, REML = TRUE), "REML")
  expect_error(frm(bfm, dd, importance = 100L,
                   control = frmtmb_control(profile = TRUE)),
               "profile")
  expect_error(frm(bfm, dd, importance = 100L, dry_run = "objective"),
               "no unfitted form")
  expect_error(frmtmb_control(importance_ess = 2), "FRACTION")
  expect_error(frmtmb_control(importance_rounds = 0), "importance_rounds")
  expect_error(frmtmb_control(importance_seed = -1), "importance_seed")
})

test_that("importance refuses the model structures it cannot correct", {
  dd <- imp_gauss_data()
  # no random-effect block at all
  expect_error(frm(bf(y ~ x) + gaussian(), dd, importance = 100L),
               "exactly one random-effect block")
  # two blocks make the integral a nested one
  dd2 <- dd
  dd2$h <- factor(rep(seq_len(10), length.out = nrow(dd2)))
  expect_error(frm(bf(y ~ x + (1 | g) + (1 | h)) + gaussian(), dd2,
                   importance = 100L),
               "exactly one random-effect block")
  # a smooth is one field over all rows, not a set of groups
  expect_error(frm(bf(y ~ s(x)) + gaussian(), dd, importance = 100L),
               "grouping factor")
  # a residual correlation term is a joint density per group
  set.seed(7)
  ad <- data.frame(y = rnorm(80), x = rnorm(80),
                   wk = rep(1:8, 10),
                   subj = factor(rep(1:10, each = 8)),
                   g = factor(rep(1:10, each = 8)))
  expect_error(
    frm(bf(y ~ x + (1 | g) + ar(wk, subj, cov = TRUE)) + gaussian(), ad,
        importance = 100L),
    "residual correlation")
})

test_that("a multi-membership term refuses: its rows span groups", {
  set.seed(31)
  n <- 200L
  dd <- data.frame(x = rnorm(n),
                   g1 = factor(sample(1:12, n, TRUE)),
                   g2 = factor(sample(1:12, n, TRUE)))
  dd$y <- rnorm(n, 0.5 * dd$x, 1)
  expect_error(frm(bf(y ~ x + (1 | mm(g1, g2))) + gaussian(), dd,
                   importance = 100L),
               "one grouping level")
})

test_that("a nonlinear predictor refuses before any tape is built", {
  # Checked through the scope guard rather than through frm(), so the
  # test does not also depend on finding starting values a nonlinear
  # body is defined at.
  set.seed(41)
  nd <- data.frame(t = rep(seq(0, 10, length.out = 25), 4),
                   g = factor(rep(1:4, each = 25)))
  nd$y <- 8 * (1 - exp(-nd$t / 3)) + rnorm(100, 0, 0.3)
  bfm <- bf(y ~ asym * (1 - exp(-t / lrc)),
            asym ~ 1 + (1 | g), lrc ~ 1, nl = TRUE) + gaussian()
  fr <- frm(bfm, data = nd, dry_run = "frame")
  tpl <- frmtmb:::make_start(fr, NULL, NULL)
  expect_error(
    frmtmb:::check_importance_scope(fr[["spec"]], fr, tpl, FALSE, FALSE,
                                    frmtmb_control()),
    "nonlinear predictor")
})

# ---------------------------------------------------------------------
# (h) DETERMINISM, and RNG hygiene.
# ---------------------------------------------------------------------
test_that("the same seed gives the same answer, and the session RNG is safe", {
  dd <- imp_scalar_data()
  bfm <- bf(y ~ x + (1 | g)) + binomial()
  f1 <- frm(bfm, dd, importance = 200L)
  f2 <- frm(bfm, dd, importance = 200L)
  expect_identical(f1$opt$objective, f2$opt$objective)
  expect_identical(f1$opt$par, f2$opt$par)
  # a different seed gives a different answer, but only by Monte Carlo
  f3 <- frm(bfm, dd, importance = 200L,
            control = frmtmb_control(importance_seed = 99L))
  expect_false(identical(f1$opt$objective, f3$opt$objective))
  expect_lt(abs(f1$opt$objective - f3$opt$objective),
            5 * f1$importance$mcse)
  # the draws come from a private stream: a fit neither reads nor
  # disturbs the session's random state
  set.seed(99)
  before <- rnorm(3)
  invisible(frm(bfm, dd, importance = 200L))
  set.seed(99)
  expect_identical(before, rnorm(3))
})

# ---------------------------------------------------------------------
# What the correction supports that quadrature has to refuse, because
# its draws sit near the conditional mode instead of out in the tails.
# ---------------------------------------------------------------------
test_that("trunc() and cens() are corrected, not refused", {
  set.seed(21)
  ng <- 40L
  per <- 6L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u <- rnorm(ng, 0, 0.8)
  y <- rpois(ng * per, exp(1.1 + 0.4 * x + u[g]))
  dd <- data.frame(y = y, x = x, g = g)[y >= 1L, ]
  bfm <- bf(y | trunc(lb = 1) ~ x + (1 | g)) + poisson()
  # quadrature cannot: its far nodes underflow the normalizer
  expect_error(frm(bfm, dd, quadrature = TRUE), "trunc")
  fi <- frm(bfm, dd, importance = 500L)
  expect_true(is.finite(fi$opt$objective))
  expect_gt(fi$importance$ess_min, frmtmb:::imp_ess_floor)

  set.seed(22)
  tt <- exp(1 + 0.5 * x + u[g] + rnorm(ng * per, 0, 0.5))
  cd <- data.frame(y = pmin(tt, 6), cen = ifelse(tt > 6, 1, 0),
                   x = x, g = g)
  fc <- frm(bf(y | cens(cen) ~ x + (1 | g)) + lognormal(), cd,
            importance = 500L)
  expect_true(is.finite(fc$opt$objective))
  expect_gt(fc$importance$ess_min, frmtmb:::imp_ess_floor)
})

# ---------------------------------------------------------------------
# What the fit records, and what it says about itself.
# ---------------------------------------------------------------------
test_that("the fit records the correction, and print/summary say so", {
  dd <- imp_scalar_data()
  fi <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd,
            importance = 200L)
  im <- fi$importance
  expect_identical(im$draws, 200L)
  expect_identical(im$seed, 1L)
  expect_gte(im$rounds, 1L)
  expect_true(is.finite(im$mcse))
  expect_length(im$ess, 100L)
  expect_identical(names(im$ess), levels(dd$g))
  # an odd request is rounded UP so the antithetic draws pair
  fo <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd,
            importance = 51L)
  expect_identical(fo$importance$draws, 52L)
  out <- paste(utils::capture.output(print(fi)), collapse = " ")
  expect_match(out, "importance-corrected")
  sout <- paste(utils::capture.output(print(summary(fi))), collapse = " ")
  expect_match(sout, "importance-corrected")
  # a Laplace fit says nothing about a correction it did not make
  fl <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd)
  expect_null(fl$importance)
  lout <- paste(utils::capture.output(print(fl)), collapse = " ")
  expect_false(grepl("importance-corrected", lout))
})

test_that("the per-group pieces reproduce the joint density they sum", {
  # imp_verify() is the pin that keeps the stacked reimplementation from
  # drifting away from the objective it corrects.
  p <- imp_parts(bf(y ~ x + (1 | g)) + binomial(), imp_scalar_data())
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$bk, p$opt$par, 50L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$bk, p$gmap, plan)
  expect_silent(frmtmb:::imp_verify(io, p$nll, plan, p$tpl, p$opt$par))
  # every row is counted exactly once, in exactly one group
  expect_length(p$gmap[["row_level"]], 400L)
  expect_identical(sort(unique(p$gmap[["row_level"]])), seq_len(100L))
  expect_equal(as.numeric(Matrix::rowSums(p$gmap[["S"]])), rep(4, 100),
               tolerance = 1e-12)
})

# An unidentified covariance turns the fixed-point iteration into a
# wander: the rounds chase each other out to a collapsed variance
# component where every weight is trivially equal, so the effective
# sample sizes at the END look perfect while the estimate is nonsense.
# The invariant that catches it is that the correction must not make
# the fit worse than the correction at the Laplace estimates.
test_that("a diverging iteration refuses instead of reporting nonsense", {
  set.seed(11)
  ng <- 40L
  per <- 4L
  g <- rep(seq_len(ng), each = per)
  x <- rnorm(ng * per)
  b <- cbind(rnorm(ng, 0, 1), rnorm(ng, 0, 0.7))
  dd <- data.frame(y = rbinom(ng * per, 1,
                              plogis(-0.5 + 0.8 * x + b[g, 1] +
                                       b[g, 2] * x)),
                   x = x, g = factor(g))
  # four binary rows per group cannot identify a 2 x 2 covariance, and
  # the Laplace fit says so itself
  fl <- suppressWarnings(frm(bf(y ~ x + (x | g)) + binomial(), data = dd))
  expect_gt(fl$opt$convergence, 0)
  expect_error(
    suppressWarnings(frm(bf(y ~ x + (x | g)) + binomial(), data = dd,
                         importance = 500L)),
    "ROSE from")
  # and the correction AT the Laplace estimates is still sound: it is
  # the iteration that fails, not the estimator
  skip_if_not_installed("GLMMadaptive")
  ga <- suppressWarnings(
    GLMMadaptive::mixed_model(y ~ x, ~ x | g, dd, binomial(), nAGQ = 25))
  p <- imp_parts(bf(y ~ x + (x | g)) + binomial(), dd)
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$bk, p$opt$par, 500L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$bk, p$gmap, plan)
  pl <- frmtmb:::imp_par_list(p$tpl, p$opt$par)
  at_lap <- io$fn(pl)
  expect_lt(abs(at_lap - -as.numeric(logLik(ga))), 0.1)
  # and it is closer to the truth than the Laplace value it corrects
  expect_lt(abs(at_lap - -as.numeric(logLik(ga))),
            abs(p$opt$objective - -as.numeric(logLik(ga))))
})

test_that("the rise tolerance separates drift from divergence", {
  # the numbers the threshold was placed from
  expect_equal(frmtmb:::imp_worse_tol(0, 0), 0.1)
  expect_equal(frmtmb:::imp_worse_tol(0.2, 0.1), 0.6)
  # the worst legitimate rise measured (gaussian optimizer drift) is
  # inside it; the smallest divergence measured is outside
  expect_lt(0.010, frmtmb:::imp_worse_tol(0, 0))
  expect_gt(0.352, frmtmb:::imp_worse_tol(0, 0))
})

# A refit must not silently drop the correction. Every path that
# rebuilds a fit already forwards `quadrature`; `importance` is
# forwarded the same way, so a profile interval or an influence refit
# optimizes the objective the fit was made with and not a different
# one.
test_that("the refit paths carry the correction rather than dropping it", {
  set.seed(601)
  ng <- 60L
  per <- 4L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u <- rnorm(ng, 0, 1.2)
  dd <- data.frame(y = rbinom(ng * per, 1, plogis(-0.5 + 0.7 * x + u[g])),
                   x = x, g = g)
  fi <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd,
            importance = 100L)
  r <- refit(fi, newresp = dd$y)
  expect_false(is.null(r$importance))
  expect_identical(r$importance$draws, fi$importance$draws)
  # same data, same seed, same objective
  expect_equal(r$opt$objective, fi$opt$objective, tolerance = 1e-6)
  # the Wald surface comes from the corrected sdreport and works
  ci <- confint(fi)
  expect_false(anyNA(ci[1:2, ]))
})

test_that("cluster-robust covariance refuses a corrected fit", {
  set.seed(601)
  ng <- 40L
  per <- 4L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u <- rnorm(ng, 0, 1.2)
  dd <- data.frame(y = rbinom(ng * per, 1, plogis(-0.5 + 0.7 * x + u[g])),
                   x = x, g = g)
  fi <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd,
            importance = 100L)
  # the corrected objective is a sum over GROUPS of reweighted
  # integrals, not the per-row sum the cluster scores are read off
  expect_error(vcov_cluster(fi, dd$g), "importance correction")
})

# P1: the compatibility table must not contradict the scope check.
# frm_compat() is what the vignette sends users to, so every structure
# the objective refuses by name has to read "refused" there, and the
# three groups together have to exhaust the registry so a structure
# added later cannot quietly report "untested".
test_that("every covariance structure importance refuses reads refused", {
  keys <- names(frmtmb:::covstruct_registry)
  ok <- frmtmb:::imp_covstructs
  cross <- frmtmb:::imp_crosslevel
  refused <- setdiff(keys, c(ok, cross))
  # the groups partition the registry: nothing unclassified
  expect_setequal(c(ok, cross, refused), keys)
  expect_length(intersect(ok, cross), 0L)
  # and the table agrees with the code, structure by structure
  for (k in c(refused, cross)) {
    expect_identical(frm_compat("importance", k)$status, "refused",
                     info = k)
  }
  for (k in ok) {
    expect_identical(frm_compat("importance", k)$status, "conditional",
                     info = k)
  }
  # the Student-t latents in particular: they factorize over levels but
  # are not gaussian in a level's coefficients, so the per-level
  # precision extraction would silently return a gaussian answer
  expect_true(all(c("us_t", "diag_t") %in% refused))
})

# P2: confint(method = "profile") does not refit an importance fit, it
# profiles the FROZEN tape, so a bound far from the estimate can come
# from a region the proposal no longer covers. It must not do that
# silently.
test_that("a profile bound in the degenerate regime warns", {
  set.seed(11)
  ng <- 40L
  per <- 4L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u <- rnorm(ng, 0, 1.1)
  dd <- data.frame(y = rbinom(ng * per, 1,
                              plogis(-0.4 + 0.7 * x + u[g])),
                   x = x, g = g)
  fi <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd,
            importance = 500L)
  # the frozen proposal, rebuilt from the fit, reproduces the anchor
  # the fit reported: same draws, same seed, same Gaussian
  prop <- frmtmb:::imp_frozen_proposal(fi)
  expect_equal(min(frmtmb:::imp_ess_at(prop, fi$opt$par)$ess),
               fi$importance$ess_min, tolerance = 1e-10)
  # the covariance parameter walks far enough to leave the proposal
  expect_warning(confint(fi, parm = "theta_1", method = "profile"),
                 "stopped covering the integrand")
  # the fixed effect does not, on this design
  expect_silent(confint(fi, parm = "x", method = "profile"))
  # and a Laplace fit profiles with none of this machinery
  fl <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd)
  expect_silent(confint(fl, parm = "theta_1", method = "profile"))
})

# P5: the pin has to be PER GROUP. A pair of errors that cancels in the
# total is exactly the kind a total-only check waves through, and the
# per-group values are what the log-mean-exp and the effective sample
# sizes rest on.
test_that("the pin catches a per-group error that cancels in the total", {
  p <- imp_parts(bf(y ~ x + (1 | g)) + binomial(), imp_scalar_data())
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$bk, p$opt$par, 50L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$bk, p$gmap, plan)
  # the honest objective passes
  expect_silent(frmtmb:::imp_verify(io, p$nll, plan, p$tpl, p$opt$par))
  # move one group up and another down by the same amount: the TOTAL is
  # untouched, so a check on the sum alone sees nothing
  bent <- io
  # in ONE draw column, so the totals at both checked columns are
  # untouched and only the per-group comparison can see it
  bent$amat <- function(pars) {
    a <- as.matrix(io$amat(pars))
    a[1L, 1L] <- a[1L, 1L] + 5
    a[2L, 1L] <- a[2L, 1L] - 5
    a
  }
  expect_equal(sum(as.matrix(bent$amat(
                 frmtmb:::imp_par_list(p$tpl, p$opt$par)))),
               sum(as.matrix(io$amat(
                 frmtmb:::imp_par_list(p$tpl, p$opt$par)))),
               tolerance = 1e-8)
  expect_error(frmtmb:::imp_verify(bent, p$nll, plan, p$tpl, p$opt$par),
               "group")
})
