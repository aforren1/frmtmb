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
  lay <- frmtmb:::imp_layout(fr[["re_blocks"]])
  list(frame = fr, nll = nll, tpl = tpl, lap = lap, opt = opt, lay = lay,
       gmap = frmtmb:::imp_group_map(fr, lay))
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
    plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, p$opt$par, nd, 1L)
    io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
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
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, p$opt$par, 200L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
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
    plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, anchor, 1000L, 1L)
    io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
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
    frmtmb:::imp_ess_warning(at_far$ess, p$lay, 1000L,
                             frmtmb:::imp_ess_floor),
    "covers"
  )
  # by the GROUPING FACTOR and not by a term label: one proposal covers
  # a level's coefficients from every block over that factor, so there
  # is one effective sample size per level. This design has ONE block,
  # which is the case the rename also changed, so the pin belongs here.
  msg <- tryCatch(frmtmb:::imp_ess_warning(at_far$ess, p$lay, 1000L,
                                           frmtmb:::imp_ess_floor),
                  warning = conditionMessage)
  expect_match(msg, "groups of `g` poorly", fixed = TRUE)
  expect_false(grepl("x | g", msg, fixed = TRUE))
  # and the levels it lists are the factor's own labels
  expect_match(msg, "'[0-9]+' [(]")
  expect_silent(
    frmtmb:::imp_ess_warning(at_opt$ess, p$lay, 1000L,
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
  # no random-effect block at all: nothing to correct
  expect_error(frm(bf(y ~ x) + gaussian(), dd, importance = 100L),
               "needs a random-effect block to correct")
  # two blocks over DIFFERENT factors make the integral a nested one
  dd2 <- dd
  dd2$h <- factor(rep(seq_len(10), length.out = nrow(dd2)))
  expect_error(frm(bf(y ~ x + (1 | g) + (1 | h)) + gaussian(), dd2,
                   importance = 100L),
               "share ONE grouping factor")
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
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, p$opt$par, 50L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
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
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, p$opt$par, 500L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
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
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, p$opt$par, 50L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
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

# =====================================================================
# SEVERAL BLOCKS OVER ONE GROUPING FACTOR.
#
# Distributional regression writes them by construction: `(1 | g)` in
# mu and `(1 | g)` in sigma are two blocks, never one, because only an
# |ID| key merges terms. The integral still factorizes over the levels
# of g, so a level's coefficients from EVERY block are drawn together
# from one joint proposal.
#
# What that breaks is contiguity. `b` is level-major within a block but
# the blocks are concatenated, so a group owns one run per block and
# the runs sit `n_levels * dim` apart. Every test below is built to
# fail if any consumer of that layout reverts to a stride: the identity
# in (a) is algebraic and would break outright, and the reference in
# (b) shares nothing with the package at all.
# =====================================================================

# Two scalar blocks over one factor, both in mu: an uncorrelated
# intercept and slope written as two terms rather than one `(x | g)`.
imp_split_data <- function() {
  set.seed(23)
  ng <- 20L
  per <- 6L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u1 <- rnorm(ng, 0, 0.9)
  u2 <- rnorm(ng, 0, 0.6)
  data.frame(y = rnorm(ng * per, 1 + 0.5 * x + u1[g] + u2[g] * x, 1),
             x = x, g = g)
}

# `(1 | g)` in mu and `(1 | g)` in sigma: the design the brute-force
# reference below integrates by hand. G is small so a 2-D quadrature
# over every group is affordable in a test.
imp_dpar_data <- function() {
  set.seed(19)
  ng <- 12L
  per <- 10L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  b1 <- rnorm(ng, 0, 0.8)
  b2 <- rnorm(ng, 0, 0.4)
  data.frame(y = rnorm(ng * per, 0.4 + 0.6 * x + b1[g],
                       exp(-0.2 + b2[g])), x = x, g = g)
}

# ---------------------------------------------------------------------
# (a) THE GAUSSIAN IDENTITY, now with the coefficients scattered.
#
# Both blocks enter mu linearly, so for a gaussian response the joint
# negative log-density is still exactly quadratic in the whole vector
# of random effects and the Laplace Gaussian IS the conditional. Every
# weight is equal and the corrected objective returns the Laplace value
# algebraically, at any draw count.
#
# This is the sharpest test of the scattered index there is: it is an
# algebraic cancellation, so a draw placed on the wrong row, or a
# half-norm read off the wrong row, destroys it outright rather than
# degrading it.
# ---------------------------------------------------------------------
test_that("two blocks in mu: the gaussian identity survives scattering", {
  dd <- imp_split_data()
  p <- imp_parts(bf(y ~ x + (1 | g) + (0 + x | g)) + gaussian(), dd)
  ng <- nlevels(dd$g)
  # the layout really is scattered, so the identity below is testing
  # what it claims to test and not a contiguous special case
  expect_length(p$lay$blocks, 2L)
  expect_identical(p$lay$qt, 2L)
  expect_identical(p$lay$idx[, 1L], c(1L, ng + 1L))
  expect_identical(p$lay$idx[, ng], c(ng, 2L * ng))
  expect_identical(p$lay$pos_group[c(1L, ng + 1L)], c(1L, 1L))

  for (nd in c(2L, 10L, 1000L)) {
    plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, p$opt$par, nd, 1L)
    io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
    pl <- frmtmb:::imp_par_list(p$tpl, p$opt$par)
    expect_equal(io$fn(pl), p$opt$objective, tolerance = 1e-10)
    ess <- frmtmb:::imp_ess(as.matrix(io$amat(pl)), plan[["n_draw"]])
    expect_equal(min(ess$ess), 1, tolerance = 1e-8)
    expect_equal(ess$mcse, 0, tolerance = 1e-8)
  }
})

# ---------------------------------------------------------------------
# (b) A BRUTE-FORCE REFERENCE that shares nothing with the package.
#
# `(1 | g)` in mu and `(1 | g)` in sigma. Each group's marginal
# likelihood is a two-dimensional integral over (b_mu, b_sigma) of the
# gaussian density times two independent normal priors, and it is done
# here by tensor Gauss-Hermite quadrature on nodes built from scratch,
# against `dnorm()` directly. No frmtmb code is involved in the
# reference beyond reading the parameter values off the fit.
#
# Unlike the gaussian designs above, the Laplace approximation is
# genuinely WRONG here: sigma depends on its random effect through a
# log link, so the joint density is not quadratic in the random
# effects and there is a real error for the correction to remove.
# ---------------------------------------------------------------------

# Gauss-Hermite nodes and weights by Golub-Welsch: the eigenvalues of
# the Jacobi matrix of the Hermite recurrence. The squared first
# components of its eigenvectors are the weights normalized to sum to
# 1, which is exactly the normal expectation weight.
imp_gh <- function(n) {
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- sqrt(i / 2)
  J[cbind(i + 1L, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  list(x = rev(e$values), lw = log(rev(e$vectors[1L, ]^2)))
}

imp_lse <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

# -sum_g log integral over (b1, b2) of prod_i dnorm(y_i, Xb + b1,
# exp(Xs bd + b2)) with b1 ~ N(0, s1^2) and b2 ~ N(0, s2^2).
#
# nq = 150 is converged: against nq = 300 it moves by 1.1e-06, and
# against nq = 100 by 1.7e-05, both far below the 0.04 Monte Carlo
# standard error the comparison is made at.
imp_ghq_ref <- function(y, X, Xs, g, beta, betad, s1, s2, nq = 150L) {
  q <- imp_gh(nq)
  b1 <- sqrt(2) * s1 * q$x
  b2 <- sqrt(2) * s2 * q$x
  lwo <- outer(q$lw, q$lw, "+")
  eta <- as.numeric(X %*% beta)
  etas <- as.numeric(Xs %*% betad)
  tot <- 0
  for (k in levels(g)) {
    ii <- which(g == k)
    m <- length(ii)
    # sigma varies with the b_sigma node only, so one m x nq matrix
    # serves every b_mu node
    sig <- exp(outer(etas[ii], b2, "+"))
    lg <- matrix(0, nq, nq)
    for (a in seq_len(nq)) {
      lg[a, ] <- colSums(matrix(stats::dnorm(y[ii], eta[ii] + b1[a], sig,
                                             log = TRUE), m, nq))
    }
    tot <- tot + imp_lse(as.numeric(lg + lwo))
  }
  -tot
}

test_that("mu and sigma blocks: brute force agrees, Laplace does not", {
  dd <- imp_dpar_data()
  p <- imp_parts(bf(y ~ x + (1 | g), sigma ~ 1 + (1 | g)) + gaussian(), dd)
  expect_length(p$lay$blocks, 2L)
  # group 1 owns one position in each block, and they are not adjacent
  expect_identical(p$lay$idx[, 1L], c(1L, nlevels(dd$g) + 1L))

  pl <- frmtmb:::imp_par_list(p$tpl, p$opt$par)
  ref <- imp_ghq_ref(dd$y, model.matrix(~ x, dd), model.matrix(~ 1, dd),
                     dd$g, pl$beta, pl$betad, exp(pl$theta[1L]),
                     exp(pl$theta[2L]))

  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, p$opt$par, 2000L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
  ess <- frmtmb:::imp_ess(as.matrix(io$amat(pl)), plan[["n_draw"]])
  got <- io$fn(pl)

  # measured: reference 193.161829, importance 193.137101 (0.61 MCSE
  # away, MCSE 0.0404), Laplace 193.396458 (5.8 MCSE away). Over ten
  # seeds the worst importance miss was 1.30 MCSE, so 3 has a 2.3x
  # margin over what was actually observed.
  expect_lt(abs(got - ref), 3 * ess$mcse)
  # and the Laplace error is real, not a rounding difference: this is
  # the whole reason the correction exists
  expect_gt(abs(p$opt$objective - ref), 4 * ess$mcse)
  # the proposal covers this design comfortably at every level
  expect_gt(min(ess$ess), frmtmb:::imp_ess_floor)
  # the per-group pieces still reproduce the joint density
  expect_silent(frmtmb:::imp_verify(io, p$nll, plan, p$tpl, p$opt$par))
})

test_that("a fit over mu and sigma blocks records its correction", {
  dd <- imp_dpar_data()
  bfm <- bf(y ~ x + (1 | g), sigma ~ 1 + (1 | g)) + gaussian()
  fl <- frm(bfm, data = dd)
  fi <- frm(bfm, data = dd, importance = 500L)
  im <- fi$importance
  expect_identical(im$draws, 500L)
  expect_false(im$capped)
  # one effective sample size per LEVEL, not per block
  expect_length(im$ess, nlevels(dd$g))
  expect_identical(names(im$ess), levels(dd$g))
  expect_gt(im$ess_min, frmtmb:::imp_ess_floor)
  # the correction moved the answer, and in the direction it must:
  # the Laplace objective here is biased upward
  expect_lt(-as.numeric(logLik(fi)), -as.numeric(logLik(fl)))
  out <- paste(utils::capture.output(print(fi)), collapse = " ")
  expect_match(out, "importance-corrected")
})

# ---------------------------------------------------------------------
# (b2) THE SAME ARBITER, WITHOUT THE GAUSSIAN.
#
# Every other multi-block check here is gaussian somewhere: (a) is an
# algebraic identity that holds because the joint is quadratic, and (b)
# integrates a gaussian response. The one non-gaussian multi-block test
# is a numDeriv gradient comparison, which pins the DERIVATIVE of the
# objective and not its value.
#
# So this is the value check with nothing gaussian left in it: 30
# groups of 8 Bernoulli rows over two blocks, integrated by the same
# tensor Gauss-Hermite construction against `dbinom()` directly. The
# quantity is fixed by no identity, and the Laplace approximation has a
# real error to remove, so a draw placed on the wrong group's rows
# cannot cancel out of it.
#
# Runs in about 6 seconds, under the file-level skip_on_cran().
# ---------------------------------------------------------------------

# Two scalar blocks, `(1 | g)` and `(0 + x | g)`, so group 1 owns
# positions 1 and ng + 1. Binary rows in small clusters are where the
# Laplace approximation is worst, which is what makes the comparison
# below worth making.
imp_bern_data <- function() {
  set.seed(5)
  ng <- 30L
  per <- 8L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u1 <- rnorm(ng, 0, 1.2)
  u2 <- rnorm(ng, 0, 0.9)
  data.frame(y = rbinom(ng * per, 1,
                        plogis(-0.4 + 0.8 * x + u1[g] + u2[g] * x)),
             x = x, g = g)
}

# -sum_g log integral over (u1, u2) of prod_i dbinom(y_i, 1,
# plogis(X beta + u1 + u2 x_i)) with u1 ~ N(0, s1^2), u2 ~ N(0, s2^2).
# Nodes from `imp_gh()` above; no frmtmb code beyond reading the
# parameter values off the fit.
#
# nq = 70 is converged: against nq = 100 it moves by 1e-11 and against
# nq = 50 by 8e-09, both far below the Monte Carlo standard error the
# comparison is made at.
imp_ghq_bern_ref <- function(y, X, x, g, beta, s1, s2, nq = 70L) {
  q <- imp_gh(nq)
  u1 <- sqrt(2) * s1 * q$x
  u2 <- sqrt(2) * s2 * q$x
  lwo <- outer(q$lw, q$lw, "+")
  eta <- as.numeric(X %*% beta)
  tot <- 0
  for (k in levels(g)) {
    ii <- which(g == k)
    m <- length(ii)
    # the slope node enters through x, so one m x nq matrix of slope
    # contributions serves every intercept node
    sl <- outer(x[ii], u2)
    lg <- matrix(0, nq, nq)
    for (a in seq_len(nq)) {
      pr <- stats::plogis(eta[ii] + u1[a] + sl)
      lg[a, ] <- colSums(matrix(stats::dbinom(y[ii], 1, pr, log = TRUE),
                                m, nq))
    }
    tot <- tot + imp_lse(as.numeric(lg + lwo))
  }
  -tot
}

test_that("two blocks, Bernoulli: brute force agrees, Laplace does not", {
  dd <- imp_bern_data()
  p <- imp_parts(bf(y ~ x + (1 | g) + (0 + x | g)) + binomial(), dd)
  ng <- nlevels(dd$g)
  expect_length(p$lay$blocks, 2L)
  # the coefficients really are scattered on this design too
  expect_identical(p$lay$idx[, 1L], c(1L, ng + 1L))

  pl <- frmtmb:::imp_par_list(p$tpl, p$opt$par)
  X <- model.matrix(~ x, dd)
  ref <- imp_ghq_bern_ref(dd$y, X, dd$x, dd$g, pl$beta,
                          exp(pl$theta[1L]), exp(pl$theta[2L]))
  # the reference is pinned before it is trusted: a coarser rule gives
  # the same number
  expect_equal(imp_ghq_bern_ref(dd$y, X, dd$x, dd$g, pl$beta,
                                exp(pl$theta[1L]), exp(pl$theta[2L]),
                                nq = 50L),
               ref, tolerance = 1e-7)

  # measured: reference 140.180284, Laplace 141.040486 (off by 0.860),
  # importance 140.175338 at 2000 draws (0.17 MCSE away, MCSE 0.0296)
  # and 140.184317 at 8000 (0.26 MCSE, MCSE 0.0154, so the Laplace
  # error is 56 MCSE there). The correction removes 99.5 percent of an
  # error that no identity fixes.
  for (nd in c(2000L, 8000L)) {
    plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, p$opt$par, nd, 1L)
    io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
    ess <- frmtmb:::imp_ess(as.matrix(io$amat(pl)), plan[["n_draw"]])
    expect_lt(abs(io$fn(pl) - ref), 3 * ess$mcse)
    expect_gt(abs(p$opt$objective - ref), 10 * ess$mcse)
    # and the proposal covers this design, so the agreement is not a
    # lucky cancellation in a degenerate one
    expect_gt(min(ess$ess), frmtmb:::imp_ess_floor)
  }
})

# ---------------------------------------------------------------------
# (c) A DIM-2 BLOCK PLUS A SCALAR ONE, the case where a group's
# coefficients are scattered AND its Hessian slice is a genuine 3 x 3
# with cross-block entries. The gradient of the corrected objective is
# checked against numDeriv, and the proposal's coverage at the optimum
# is reported.
# ---------------------------------------------------------------------
test_that("(1 + x | g) plus a sigma block: gradient and coverage hold", {
  skip_if_not_installed("numDeriv")
  dd <- imp_dpar_data()
  p <- imp_parts(bf(y ~ x + (1 + x | g), sigma ~ 1 + (1 | g)) + gaussian(),
                 dd)
  ng <- nlevels(dd$g)
  expect_identical(p$lay$qt, 3L)
  # the mu block's two coefficients are adjacent; the sigma block's
  # sits 2 * ng away, which is the whole point
  expect_identical(p$lay$idx[, 1L], c(1L, 2L, 2L * ng + 1L))
  expect_identical(p$lay$idx[, 2L], c(3L, 4L, 2L * ng + 2L))

  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, p$opt$par, 500L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
  expect_silent(frmtmb:::imp_verify(io, p$nll, plan, p$tpl, p$opt$par))

  otpl <- frmtmb:::imp_template(p$tpl, "b")
  obj <- RTMB::MakeADFun(io$fn, otpl, silent = TRUE)
  # The tolerance is set by numDeriv, not by the AD gradient: two
  # numDeriv settings disagree with each other by 3.0e-04 at the
  # optimum and 1.0e-02 away from it, which is larger than the gap
  # being measured (relative 1.6e-06 and 6.9e-06).
  for (shift in c(0, 0.3)) {
    at <- p$opt$par + shift
    nd <- numDeriv::grad(obj$fn, at)
    expect_equal(as.numeric(obj$gr(at)), nd, tolerance = 1e-4)
  }
  ess <- frmtmb:::imp_ess(
    as.matrix(io$amat(frmtmb:::imp_par_list(p$tpl, p$opt$par))),
    plan[["n_draw"]])
  # measured: min 0.172, median 0.881 at 500 draws
  expect_gt(min(ess$ess), 0.05)
  expect_gt(stats::median(ess$ess), frmtmb:::imp_ess_floor)
})

# ---------------------------------------------------------------------
# (d) THE REFUSALS several blocks add.
# ---------------------------------------------------------------------
test_that("blocks over different factors are refused by name", {
  dd <- imp_gauss_data()
  dd$h <- factor(rep(seq_len(10), length.out = nrow(dd)))
  err <- tryCatch(frm(bf(y ~ x + (1 | g) + (1 | h)) + gaussian(), dd,
                      importance = 100L),
                  error = conditionMessage)
  expect_match(err, "share ONE grouping factor")
  # the message names the blocks AND the factors, so the user can see
  # which term is the problem
  expect_match(err, "`1 \\| g` over g")
  expect_match(err, "`1 \\| h` over h")
  # nesting is the same refusal: still two factors
  dd$n <- factor(paste0(dd$g, ".", rep(1:2, length.out = nrow(dd))))
  expect_error(frm(bf(y ~ x + (1 | g) + (1 | n)) + gaussian(), dd,
                   importance = 100L),
               "share ONE grouping factor")
})

test_that("blocks over one factor with different level sets are refused", {
  # No user-facing route builds this: two blocks over the same column
  # take their levels from that column and so always agree. The guard
  # is checked directly on a doctored frame, because a level set that
  # silently disagreed would put a group's draws on another group's
  # rows rather than raise anything.
  dd <- imp_dpar_data()
  fr <- frm(bf(y ~ x + (1 | g), sigma ~ 1 + (1 | g)) + gaussian(),
            data = dd, dry_run = "frame")
  tpl <- frmtmb:::make_start(fr, NULL, NULL)
  expect_silent(frmtmb:::check_importance_scope(fr$spec, fr, tpl, FALSE,
                                                FALSE, frmtmb_control()))
  bad <- fr
  bad$re_blocks[[2L]]$levels <- bad$re_blocks[[2L]]$levels[-3L]
  err <- tryCatch(frmtmb:::check_importance_scope(bad$spec, bad, tpl, FALSE,
                                                  FALSE, frmtmb_control()),
                  error = conditionMessage)
  expect_match(err, "same grouping levels")
  expect_match(err, "first difference: '3'")
})

# ---------------------------------------------------------------------
# (e) THE PER-GROUP PIN still fires when a group is scattered. A bug
# that put one block's draw on the wrong group's rows would show up
# here, and this is the check that says so.
# ---------------------------------------------------------------------
test_that("the per-group pin fires with several blocks", {
  dd <- imp_dpar_data()
  p <- imp_parts(bf(y ~ x + (1 | g), sigma ~ 1 + (1 | g)) + gaussian(), dd)
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, p$opt$par, 50L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
  expect_silent(frmtmb:::imp_verify(io, p$nll, plan, p$tpl, p$opt$par))
  # every row belongs to exactly one group, read across BOTH blocks
  expect_identical(sort(unique(p$gmap[["row_level"]])),
                   seq_len(nlevels(dd$g)))
  expect_equal(as.numeric(Matrix::rowSums(p$gmap[["S"]])),
               rep(10, nlevels(dd$g)), tolerance = 1e-12)
  # and an error that cancels in the total is still caught per group
  bent <- io
  bent$amat <- function(pars) {
    a <- as.matrix(io$amat(pars))
    a[1L, 1L] <- a[1L, 1L] + 5
    a[2L, 1L] <- a[2L, 1L] - 5
    a
  }
  expect_error(frmtmb:::imp_verify(bent, p$nll, plan, p$tpl, p$opt$par),
               "group")
})

# ---------------------------------------------------------------------
# (f) DETERMINISM, with the draws spread over two blocks.
# ---------------------------------------------------------------------
test_that("several blocks stay deterministic at one seed", {
  dd <- imp_dpar_data()
  bfm <- bf(y ~ x + (1 | g), sigma ~ 1 + (1 | g)) + gaussian()
  # 200 draws hits the round cap on this design and warns about it.
  # That is convergence, not determinism: what is under test is that
  # the draw stream is a function of the seed alone.
  f1 <- suppressWarnings(frm(bfm, data = dd, importance = 200L))
  f2 <- suppressWarnings(frm(bfm, data = dd, importance = 200L))
  expect_identical(f1$opt$par, f2$opt$par)
  expect_identical(f1$opt$objective, f2$opt$objective)
  expect_identical(f1$importance$ess, f2$importance$ess)
  # and the session's own stream is untouched by a fit that draws
  set.seed(99)
  before <- rnorm(3)
  set.seed(99)
  invisible(suppressWarnings(frm(bfm, data = dd, importance = 200L)))
  expect_identical(rnorm(3), before)
})
