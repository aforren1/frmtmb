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

  fl <- frm(bf(y ~ x + (x | g)) + bernoulli(), data = dd)
  fi <- frm(bf(y ~ x + (x | g)) + bernoulli(), data = dd,
            importance = 2000L, se = TRUE)
  mcse <- fi$importance$mcse
  expect_gt(mcse, 0)

  # the corrected log-likelihood, against its own reported MCSE
  expect_lt(abs(fi$opt$objective - ref_nll), 3 * mcse)
  # the intercept, against its own standard error
  se_int <- sqrt(diag(vcov(fi)))[[1L]]
  expect_lt(abs(fixef_by_dpar(fi)$mu[[1L]] - ref_fix[[1L]]), 0.5 * se_int)
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
  fq <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd,
            quadrature = TRUE)
  fi <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd,
            importance = 1000L)
  ref <- lme4::glmer(y ~ x + (1 | g), dd, family = binomial, nAGQ = 25)

  mcse <- fi$importance$mcse
  expect_lt(abs(as.numeric(logLik(fi)) - as.numeric(logLik(ref))),
            3 * mcse)
  expect_lt(abs(as.numeric(logLik(fi)) - as.numeric(logLik(fq))),
            3 * mcse)
  expect_lt(max(abs(fixef_by_dpar(fi)$mu - lme4::fixef(ref))), 0.05)
  sd_i <- sqrt(varcorr_matrices(fi)[[1L]][1, 1])
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
  fl <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd)
  fi <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd,
            importance = 2000L)
  sd_lap <- sqrt(varcorr_matrices(fl)[[1L]][1, 1])
  sd_imp <- sqrt(varcorr_matrices(fi)[[1L]][1, 1])
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
  p <- imp_parts(bf(y ~ x + (1 | g)) + bernoulli(), imp_scalar_data())
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
# (f) THE ESS DIAGNOSTIC fires on a proposal read deliberately far
# from the optimum and stays quiet at the optimum.
# ---------------------------------------------------------------------
test_that("the effective sample size diagnostic separates the two regimes", {
  skip_if_not_installed("MASS")
  p <- imp_parts(bf(y ~ x + (x | g)) + bernoulli(), imp_probe_data())
  ess_at <- function(anchor, at) {
    plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, anchor, 1000L, 1L)
    io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
    frmtmb:::imp_ess(as.matrix(io$amat(frmtmb:::imp_par_list(p$tpl, at))),
                     plan[["n_draw"]])
  }
  at_opt <- ess_at(p$opt$par, p$opt$par)
  expect_gt(min(at_opt$ess), frmtmb:::imp_ess_floor)
  # the same draws read where they stop covering: three log standard
  # deviations up, the worst group falls under the floor and the error the
  # fit would report rises with it
  far <- p$opt$par
  far[3:4] <- far[3:4] + 3
  at_far <- ess_at(p$opt$par, far)
  expect_lt(min(at_far$ess), frmtmb:::imp_ess_floor)
  expect_gt(at_far$mcse, 5 * at_opt$mcse)
  # a design that defeats the proposal for real: 30 groups of 2 null
  # binary rows under a 2 x 2 covariance, whose group '20' has weights
  # with a tail the Laplace Gaussian does not cover, at this data seed. The
  # estimate stays unbiased, so what degrades is efficiency, and that is
  # exactly what this diagnostic is for.
  set.seed(5)
  ng <- 30L
  g <- factor(rep(seq_len(ng), each = 2L))
  xu <- rnorm(2L * ng)
  du <- data.frame(y = rbinom(2L * ng, 1, 0.4), x = xu, g = g)
  allow_warnings(
    bad <- frm(bf(y ~ x + (x | g)) + bernoulli(), data = du,
               importance = 500L),
    c("covers", "false convergence"), require = "covers")
  ess_bad <- bad$importance$ess
  expect_lt(min(ess_bad), frmtmb:::imp_ess_floor)
  # the rest of the groups hold, so it is that group the warning names
  expect_gt(stats::median(ess_bad), frmtmb:::imp_ess_floor)
  # by the GROUPING FACTOR and not by a term label: one proposal covers
  # a level's coefficients from every block over that factor, so there
  # is one effective sample size per level. This design has ONE block,
  # which is the case the rename also changed, so the pin belongs here.
  lay_bad <- frmtmb:::imp_layout(bad$frame[["re_blocks"]])
  msg <- tryCatch(frmtmb:::imp_ess_warning(ess_bad, lay_bad, 500L,
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
  bfm <- bf(y ~ x + (1 | g)) + bernoulli()
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
  fi <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd,
            importance = 200L)
  im <- fi$importance
  expect_identical(im$draws, 200L)
  expect_identical(im$seed, 1L)
  expect_gte(im$rounds, 1L)
  expect_true(is.finite(im$mcse))
  expect_length(im$ess, 100L)
  expect_identical(names(im$ess), levels(dd$g))
  # an odd request is rounded UP so the antithetic draws pair
  fo <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd,
            importance = 51L)
  expect_identical(fo$importance$draws, 52L)
  out <- paste(utils::capture.output(print(fi)), collapse = " ")
  expect_match(out, "importance-corrected")
  sout <- paste(utils::capture.output(print(summary(fi))), collapse = " ")
  expect_match(sout, "importance-corrected")
  # a Laplace fit says nothing about a correction it did not make
  fl <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd)
  expect_null(fl$importance)
  lout <- paste(utils::capture.output(print(fl)), collapse = " ")
  expect_false(grepl("importance-corrected", lout))
})

test_that("the per-group pieces reproduce the joint density they sum", {
  # imp_verify() is the pin that keeps the stacked reimplementation from
  # drifting away from the objective it corrects.
  p <- imp_parts(bf(y ~ x + (1 | g)) + bernoulli(), imp_scalar_data())
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, p$opt$par, 50L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
  expect_silent(frmtmb:::imp_verify(io, p$nll, plan, p$tpl, p$opt$par))
  # every row is counted exactly once, in exactly one group
  expect_length(p$gmap[["row_level"]], 400L)
  expect_identical(sort(unique(p$gmap[["row_level"]])), seq_len(100L))
  expect_equal(as.numeric(Matrix::rowSums(p$gmap[["S"]])), rep(4, 100),
               tolerance = 1e-12)
})
# Four binary rows per group, a 2 x 2 covariance, a Laplace fit that stops
# at the correlation's boundary with a convergence code. This design used
# to stop with "ROSE from" at 500 draws, and the message blamed an
# unidentified covariance. The cause was the frozen proposal's noisy
# curvature: with the proposal following theta, it settles on the exact
# maximum.
test_that("a weakly identified covariance settles on the exact value", {
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
  fl <- suppressWarnings(frm(bf(y ~ x + (x | g)) + bernoulli(), data = dd))
  expect_gt(fl$opt$convergence, 0)
  fi <- suppressWarnings(frm(bf(y ~ x + (x | g)) + bernoulli(), data = dd,
                             importance = 500L))
  expect_false(fi$importance$capped)
  expect_gt(as.numeric(logLik(fi)), as.numeric(logLik(fl)))
  skip_if_not_installed("GLMMadaptive")
  ga <- suppressWarnings(
    GLMMadaptive::mixed_model(y ~ x, ~ x | g, dd, binomial(), nAGQ = 25))
  expect_lt(abs(as.numeric(logLik(fi)) - as.numeric(logLik(ga))),
            3 * fi$importance$mcse)
  # the estimates too, on the scale of their own standard errors
  se_b <- sqrt(diag(vcov(fi)))
  b_fi <- fi$opt$par[names(fi$opt$par) == "beta"]
  expect_lt(max(abs(b_fi - GLMMadaptive::fixef(ga)) / se_b[seq_along(b_fi)]),
            0.1)
  # and the correction AT the Laplace estimates is sound too
  p <- imp_parts(bf(y ~ x + (x | g)) + bernoulli(), dd)
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
  # the floor, and the Monte Carlo term over it
  expect_equal(frmtmb:::imp_worse_tol(0, 0), 0.1)
  expect_equal(frmtmb:::imp_worse_tol(0.2, 0.1), 0.6)
  # what a settled fit's rise reads against it, on a binary design with room
  # to rise. Over 648 measured fits the largest was 0.15 of the tolerance.
  fit <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = imp_scalar_data(),
             importance = 200L)
  rise <- fit$opt$objective - fit$importance$start_value
  tol <- frmtmb:::imp_worse_tol(fit$importance$mcse, fit$importance$mcse)
  expect_lt(abs(rise) / tol, 1)
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
  fi <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd,
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
  fi <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd,
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
  keys <- frmtmb:::cs_names()
  ok <- frmtmb:::imp_covstructs
  cross <- frmtmb:::imp_crosslevel()
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
# profiles its own tape, so a bound far from the estimate can come
# from a region the proposal no longer covers. It must not do that
# silently.
test_that("a profile bound warns only where the weights fail", {
  set.seed(11)
  ng <- 40L
  per <- 4L
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u <- rnorm(ng, 0, 1.1)
  dd <- data.frame(y = rbinom(ng * per, 1,
                              plogis(-0.4 + 0.7 * x + u[g])),
                   x = x, g = g)
  fi <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd,
            importance = 500L)
  # the proposal, rebuilt from the fit, reproduces the anchor the fit
  # reported: same draws, same seed, same Gaussian
  prop <- frmtmb:::imp_fit_proposal(fi)
  expect_equal(min(frmtmb:::imp_ess_at(prop, fi$opt$par)$ess),
               fi$importance$ess_min, tolerance = 1e-10)
  # the covariance parameter walks to bounds where the proposal, which
  # follows it, still covers: no warning, and the bounds are quadrature's
  # within a hundredth of the interval. A frozen proposal warned here.
  ci <- expect_silent(confint(fi, parm = "theta_1", method = "profile"))
  qd <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd, quadrature = TRUE)
  cq <- confint(qd, parm = "theta_1", method = "profile")
  expect_lt(max(abs(ci[1, 1:2] - cq[1, 1:2])) / (cq[1, 2] - cq[1, 1]), 0.01)
  expect_silent(confint(fi, parm = "x", method = "profile"))
  # a bound past where it covers does warn, and it names the parameter
  i <- which(names(fi$opt$par) == "theta")
  far <- fi$opt$par[[i]] + 3
  expect_warning(frmtmb:::imp_profile_ess_warn(fi, prop, "theta_1", i, far),
                 "stopped covering the integrand")
  expect_silent(frmtmb:::imp_profile_ess_warn(fi, prop, "theta_1", i,
                                              ci[1, 1:2]))
  # and a Laplace fit profiles with none of this machinery
  fl <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd)
  expect_silent(confint(fl, parm = "theta_1", method = "profile"))
})

# P5: the pin has to be PER GROUP. A pair of errors that cancels in the
# total is exactly the kind a total-only check waves through, and the
# per-group values are what the log-mean-exp and the effective sample
# sizes rest on.
test_that("the pin catches a per-group error that cancels in the total", {
  p <- imp_parts(bf(y ~ x + (1 | g)) + bernoulli(), imp_scalar_data())
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
  p <- imp_parts(bf(y ~ x + (1 | g) + (0 + x | g)) + bernoulli(), dd)
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
  # numDeriv's defaults carry 3.0e-04 of their own error here, measured
  # against a central difference and a smaller Richardson step, which
  # both agree with the AD gradient to 3e-05
  for (shift in c(0, 0.3)) {
    at <- p$opt$par + shift
    nd <- numDeriv::grad(obj$fn, at,
                         method.args = list(eps = 1e-5, d = 1e-3))
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

# ---------------------------------------------------------------------
# (l) A STRUCTURED FAMILY THAT DECLARES HOW ITS LIKELIHOOD FACTORIZES.
#
# `loglik` returns one number for the whole response, which is why the
# correction used to refuse every structured family outright. A family
# whose likelihood does factorize over the grouping says so with
# loglik_group or loglik_row, and the correction then resamples its
# pieces like any other.
#
# The family below is the strongest available control: its likelihood
# IS the gaussian one, written out. Anything the correction computes
# for it has a right answer that gaussian() gives independently.
# ---------------------------------------------------------------------

# gaussian, spelled as a structure. `pieces` chooses which factorization
# it declares, so one family covers all three cases.
imp_struct_gauss <- function(pieces = c("row", "group", "none")) {
  pieces <- match.arg(pieces)
  terms <- function(y, dpars, aterms, block) {
    RTMB::dnorm(y, dpars[["mu"]], dpars[["sigma"]], log = TRUE)
  }
  st <- frmtmb_structure(
    # the grouping column belongs to no linear predictor, so the
    # family has to ask for it or the frame will not carry it
    frame_vars = function(fam) list(fam[["g_expr"]]),
    frame_block = function(resp, spec, av, mf, y, n) {
      list(group = mf[[deparse(resp$family[["g_expr"]])]], n = n)
    },
    loglik = function(y, dpars, aterms, weights, block, extra) {
      sum(weights * terms(y, dpars, aterms, block))
    },
    loglik_row = if (pieces == "row") {
      function(y, dpars, aterms, weights, block, extra) {
        weights * terms(y, dpars, aterms, block)
      }
    },
    loglik_group = if (pieces == "group") {
      function(y, dpars, aterms, weights, block, extra) {
        ll <- weights * terms(y, dpars, aterms, block)
        # replicate-major: the design is stacked once per draw, so the
        # grouping repeats with it
        g <- rep(as.integer(factor(block[["group"]])),
                 NROW(y) %/% block[["n"]])
        g <- g + rep((seq_len(NROW(y) %/% block[["n"]]) - 1L) *
                       nlevels(factor(block[["group"]])),
                     each = block[["n"]])
        as.vector(Matrix::sparseMatrix(i = g, j = seq_along(g),
                                       x = rep(1, length(g))) %*% ll)
      }
    },
    unit = "one group's rows",
    supports = list(conditional_effects = TRUE)
  )
  fam <- frmtmb_family(
    "struct_gauss", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    accepts_aterms = "weights",
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars[["mu"]], dpars[["sigma"]], log = TRUE)
    },
    init_dpars = list(mu = function(y, aterms) mean(y),
                      sigma = function(y, aterms) sd(y)),
    structure = st)
  fam[["g_expr"]] <- quote(g)
  fam
}

test_that("a family with no factorization is refused in groups, not rows", {
  dd <- imp_gauss_data()
  fam <- imp_struct_gauss("none")
  err <- expect_error(frm(bf(y ~ x + (1 | g)), data = dd, family = fam,
                          importance = 100L))
  msg <- conditionMessage(err)
  expect_match(msg, "one number for the whole response", fixed = TRUE)
  expect_match(msg, "loglik_group", fixed = TRUE)
  expect_match(msg, "loglik_row", fixed = TRUE)
  expect_match(msg, "importance = 0", fixed = TRUE)
  # and the same model fits without the correction, so what is refused
  # is the correction and not the family
  expect_s3_class(frm(bf(y ~ x + (1 | g)), data = dd, family = fam),
                  "frmtmb_fit")
})

test_that("a declared factorization is corrected, and gaussian() agrees", {
  dd <- imp_gauss_data()
  bfm <- bf(y ~ x + (1 | g))
  ref <- frm(bfm + gaussian(), data = dd, importance = 200L)
  for (kind in c("row", "group")) {
    fit <- frm(bfm, data = dd, family = imp_struct_gauss(kind),
               importance = 200L)
    # the same corrected objective, reached through the family's own
    # pieces instead of through row_lpdf()
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
                 tolerance = 1e-6)
    expect_equal(unlist(fixef_by_dpar(fit)), unlist(fixef_by_dpar(ref)),
                 tolerance = 1e-5)
    expect_equal(fit$importance$mcse, ref$importance$mcse, tolerance = 1e-6)
  }
})

test_that("the family's units must be the model's grouping levels", {
  dd <- imp_gauss_data()
  # the family groups on `g` and the model on `h`, which cuts across it:
  # a piece then belongs to two groups and summing is not defined
  dd$h <- factor(rep(seq_len(5), length.out = nrow(dd)))
  err <- expect_error(
    frm(bf(y ~ x + (1 | h)), data = dd, family = imp_struct_gauss("group"),
        importance = 100L))
  expect_match(conditionMessage(err), "same partition of the rows",
               fixed = TRUE)
  expect_match(conditionMessage(err), "50 unit(s) against 5 grouping",
               fixed = TRUE)
})

test_that("a factorization with no grouping in the block is refused", {
  dd <- imp_gauss_data()
  fam <- imp_struct_gauss("row")
  # the block builder is what carries `group`; without it the core has
  # nothing to check the model's grouping against
  fam[["structure"]][["frame_block"]] <- function(resp, spec, av, mf, y, n) {
    list(n = n)
  }
  expect_error(frm(bf(y ~ x + (1 | g)), data = dd, family = fam,
                   importance = 100L),
               "carries no `group`")
})

# ---------------------------------------------------------------------
# WHICH KIND OF CAPPED. `capped` covers two fits that want opposite
# advice, and until 0.55.0 both got the same sentence.
#
# One is short of rounds: its moves shrink by a factor of three to ten
# each round and one or two more rounds would land it. The other is
# stalled: it takes the SAME step every round, so its total shift is
# that step times the round count, and a larger cap buys a
# proportionally larger number rather than a better one.
#
# The step is a property of the DRAWS and not of the model. Twelve
# groups of three Bernoulli rows simulated with sd = 0, at 50 draws,
# walk at 0.90907 whatever the data seed, and eight groups at 50 draws
# walk at 0.949613, which is 1.3e-06 from what a frmtmb.learn fit of an
# entirely different family reports at eight subjects and 50 draws.
# Both differences are smaller than the spread between the rounds of
# either run.
# ---------------------------------------------------------------------

# sd = 0, so the Laplace fit collapses the component and the correction
# has nothing left to reweight. Twelve groups of three rows fits in
# about a fifth of a second.
imp_stall_data <- function(seed = 1L, ng = 12L, per = 3L) {
  set.seed(seed)
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  data.frame(y = rbinom(ng * per, 1, plogis(-0.3 + 0.7 * x)), x = x, g = g)
}

# every warning the expression raises, as text
imp_warnings <- function(expr) {
  ws <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws)
}
test_that("a collapsed component settles where Laplace put it", {
  # sd = 0 used to walk: the proposal stayed put, the step was a property
  # of the draws, and the fit reported 0.909 per round for five rounds and
  # a shift of 4.545. With the proposal following theta nothing moves.
  dd <- imp_stall_data()
  lap <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd)
  fit <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd, importance = 50L)
  im <- fit$importance
  expect_false(im$capped)
  se_b <- sqrt(diag(vcov(lap)))
  b <- names(fit$opt$par) == "beta"
  expect_lt(max(abs(fit$opt$par[b] - lap$opt$par[b]) / se_b[seq_len(sum(b))]),
            0.01)
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(lap))) /
              abs(as.numeric(logLik(lap))), 1e-6)
})
test_that("a correction that is still shrinking still asks for rounds", {
  # the same design as the multi-block tests, stopped one round short:
  # with a second round it lands, so a cap of one is genuinely short of
  # rounds and not stalled
  r <- imp_warnings(frm(bf(y ~ x + (1 | g)) + bernoulli(),
                        data = imp_scalar_data(), importance = 200L,
                        control = frmtmb_control(importance_rounds = 1L)))
  im <- r$value$importance
  expect_true(im$capped)
  hit <- grep("importance correction", r$warnings, value = TRUE)
  expect_length(hit, 1L)
  expect_match(hit, "still moving", fixed = TRUE)
  # one more round is the right advice: without the cap the same fit lands
  more <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = imp_scalar_data(),
              importance = 200L)
  expect_false(more$importance$capped)
})


# ---------------------------------------------------------------------
# THE PROPOSAL FOLLOWS THE OPTIMIZER'S THETA.
#
# With a proposal frozen in place, the curvature of the objective along
# an intercept is the group's conditional information times
# (1 - mean z^2) plus the marginal curvature, a sample statistic far
# larger than what it has to resolve. The re-anchored iteration then
# diverged wherever that sample statistic came out small, however good
# the weights were, and frm() stopped on the textbook models at 200 and
# 512 draws.
# ---------------------------------------------------------------------
# The log of one group's integral over u ~ N(0, s^2), from the log integrand
# of a scalar u. Centered on its peak, because integrate()'s default absolute
# tolerance accepts anything on an integrand of size e^-40, and it reads about
# 0.04 nats wrong two standard errors off the anchor of the textbook model.
imp_log_integral <- function(logf, s) {
  grid <- seq(-10 * s, 10 * s, length.out = 4001L)
  v <- vapply(grid, logf, 0)
  top <- max(v)
  keep <- range(grid[v > top - 50])
  val <- integrate(function(u) vapply(u, function(ui) exp(logf(ui) - top), 0),
                   keep[[1L]], keep[[2L]], subdivisions = 500L,
                   rel.tol = 1e-8, abs.tol = 0)$value
  top + log(val)
}

# 10 groups of 20, rows in group order: the models that stopped.
imp_textbook_data <- function(seed) {
  set.seed(seed)
  n <- 200L
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  g = factor(sprintf("g%02d", sort(rep(1:10, length.out = n)))))
  # the draws two more columns take from this stream in the survey that
  # found the stop, so that these data are its data
  sample(c("a", "b", "c"), n, TRUE)
  sample(0:3, n, TRUE)
  u <- rnorm(10L, 0, 0.5)
  v <- rnorm(10L, 0, 0.3)
  d$eta <- 0.3 + 0.5 * d$x - 0.3 * d$z + u[as.integer(d$g)] +
    v[as.integer(d$g)] * d$x
  d$ntr <- 10L
  d
}
test_that("the textbook binomial settles at the draw counts users pick", {
  d <- imp_textbook_data(124L)
  d$y <- rbinom(nrow(d), 10L, plogis(d$eta))
  f <- bf(y | trials(ntr) ~ x + (1 | g))
  lap <- frm(f, family = binomial(), data = d)
  qd <- frm(f, family = binomial(), data = d, quadrature = TRUE)
  se_b <- sqrt(diag(vcov(qd)))
  gap_lap <- abs(as.numeric(logLik(lap)) - as.numeric(logLik(qd)))
  for (k in c(200L, 512L, 1000L)) {
    fit <- frm(f, family = binomial(), data = d, importance = k)
    im <- fit$importance
    expect_false(im$capped)
    # one round and one to confirm at most: nothing is left to iterate on
    expect_lte(im$rounds, 2L)
    # within its own reported error of the exact value, and nearer than
    # Laplace, whose gap is 0.011
    gap <- abs(as.numeric(logLik(fit)) - as.numeric(logLik(qd)))
    expect_lt(gap, 3 * im$mcse)
    expect_lt(gap, gap_lap)
    # and the estimates, on the scale of their own standard errors
    b <- fit$opt$par[names(fit$opt$par) == "beta"]
    bq <- qd$opt$par[names(qd$opt$par) == "beta"]
    expect_lt(max(abs(b - bq) / se_b[seq_along(b)]), 0.05)
  }
})
test_that("a correlated slope at a boundary settles at the draw counts", {
  d <- imp_textbook_data(69L)
  d$y <- rpois(nrow(d), exp(d$eta))
  f <- bf(y ~ x + (1 + x | g))
  # At the correlation's limit the slope is a fixed multiple of the
  # intercept's draw, so the exact marginal is one integral per group
  exact_at <- function(fit) {
    vc <- VarCorr(fit)$g
    s1 <- vc$sd[1L, "Estimate"]
    s2 <- vc$sd[2L, "Estimate"]
    rho <- vc$cor[2L, "Estimate", 1L]
    b <- fit$opt$par[names(fit$opt$par) == "beta"]
    sum(vapply(split(seq_len(nrow(d)), d$g), function(r) {
      imp_log_integral(function(ui) {
        sum(dpois(d$y[r], exp(b[[1L]] + b[[2L]] * d$x[r] + ui +
                                rho * s2 / s1 * ui * d$x[r]), log = TRUE)) +
          dnorm(ui, 0, s1, log = TRUE)
      }, s1)
    }, 0))
  }
  lap <- frm(f, family = poisson(), data = d)
  se_b <- sqrt(diag(vcov(lap)))
  b_lap <- lap$opt$par[names(lap$opt$par) == "beta"]
  for (k in c(200L, 512L)) {
    fit <- frm(f, family = poisson(), data = d, importance = k)
    im <- fit$importance
    expect_false(im$capped)
    expect_lt(abs(as.numeric(logLik(fit)) - exact_at(fit)), 3 * im$mcse)
    # the old iteration walked the correlation's chart by thousands here,
    # taking the fixed effects with it
    b <- fit$opt$par[names(fit$opt$par) == "beta"]
    expect_lt(max(abs(b - b_lap) / se_b[seq_along(b)]), 0.1)
  }
})
test_that("the objective's curvature along an intercept is not a sample", {
  d <- imp_textbook_data(124L)
  d$y <- rbinom(nrow(d), 10L, plogis(d$eta))
  p <- imp_parts(bf(y | trials(ntr) ~ x + (1 | g)) + binomial(), d)
  at <- p$opt$par
  otpl <- frmtmb:::imp_template(p$tpl, "b")
  # Laplace is good to 1 part in 1000 here, so its curvature is the
  # reference, by central difference of its exact gradient
  h <- 1e-4 * max(1, abs(at[[1L]]))
  e <- replace(numeric(length(at)), 1L, h)
  exact <- (p$lap$gr(at + e) - p$lap$gr(at - e))[[1L]] / (2 * h)
  # seed 1 gave 5 and seed 14 gave -3.4 against an exact 34 at 200 draws
  for (s in c(1L, 14L)) {
    plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, at, 200L, s)
    io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
    obj <- RTMB::MakeADFun(io$fn, otpl, silent = TRUE)
    expect_lt(abs(obj$he(at)[1L, 1L] / exact - 1), 0.1)
  }
})
test_that("a gaussian fit stays at the exact optimum at any draw count", {
  # the fixed point of the iteration is the exact optimum because the
  # gradient at the anchor is exact for every draw, when the proposal's
  # curvature follows the dispersion
  set.seed(1003L)
  g <- factor(rep(1:30, each = 3L))
  x <- round(rnorm(90L), 3)
  a <- rnorm(30L)
  b <- 0.3 * a + sqrt(1 - 0.3^2) * rnorm(30L)
  d <- data.frame(y = 0.3 * x + 0.8 * a[g] + 0.4 * b[g] * x + 0.2 +
                    rnorm(90L, 0, 0.5), x = x, g = g)
  f <- bf(y ~ x + (1 + x | g))
  lap <- frm(f, family = gaussian(), data = d)
  se <- sqrt(diag(frmtmb:::sdr_of(lap)$cov.fixed))
  for (k in c(16L, 64L, 200L)) {
    fit <- frm(f, family = gaussian(), data = d, importance = k)
    expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(lap))) /
                abs(as.numeric(logLik(lap))), 1e-6)
    expect_lt(max(abs(fit$opt$par - lap$opt$par) / se), 1e-2)
  }
})
# The optimizer and the plan, wrapped so that a chosen draw count lands
# somewhere worse than where it started. The guard cannot be reached on
# demand any other way now that the iteration settles, and the wrapper
# moves only the rounds' optimum, so all the arithmetic below it is real.
imp_with_displaced_rounds <- function(below, code) {
  real_plan <- frmtmb:::imp_plan
  real_opt <- frmtmb:::optimize_obj
  last_n <- 0L
  testthat::with_mocked_bindings(
    code,
    imp_plan = function(lap_obj, frame, lay, par, n_draw, seed) {
      last_n <<- n_draw
      real_plan(lap_obj, frame, lay, par, n_draw, seed)
    },
    optimize_obj = function(obj, control, bounds = list(lower = -Inf,
                                                        upper = Inf),
                            par_units = NULL, verbose = 0L,
                            start_par = obj$par, tally = NULL) {
      o <- real_opt(obj, control, bounds, par_units, verbose, start_par, tally)
      if (last_n > 0L && last_n < below) o$par[[1L]] <- o$par[[1L]] + 1
      o
    },
    .package = "frmtmb")
}
test_that("too few draws are named with their ESS and the count that works", {
  d <- imp_scalar_data()
  f <- bf(y ~ x + (1 | g))
  # the guard's condition is absent: the same wrapper, moving nothing
  fit <- imp_with_displaced_rounds(0L, frm(f, family = bernoulli(), data = d,
                                           importance = 16L))
  expect_s3_class(fit, "frmtmb_fit")
  msg <- tryCatch(imp_with_displaced_rounds(
    64L, frm(f, family = bernoulli(), data = d, importance = 16L)),
    error = conditionMessage)
  expect_type(msg, "character")
  expect_match(msg, "`importance = 16` did not settle", fixed = TRUE)
  expect_match(msg, "effective sample size is 0", fixed = TRUE)
  expect_match(msg, "settled at importance = 64", fixed = TRUE)
  expect_match(msg, "Use importance = 64", fixed = TRUE)
  expect_false(grepl("or more", msg, fixed = TRUE))
  expect_false(grepl("ROSE", msg, fixed = TRUE))
})
test_that("a constant dpar's numeric path reads the model the tape reads", {
  # the free betad entries used to be recycled over the entries a map
  # fixes, so the offsets and the effective sample sizes read nu = sigma
  set.seed(31L)
  g <- factor(rep(1:15, each = 8L))
  x <- rnorm(120L)
  d <- data.frame(y = 0.3 + 0.5 * x + rnorm(15L, 0, 0.7)[g] +
                    0.5 * rt(120L, 5L), x = x, g = g)
  fit <- frm(bf(y ~ x + (1 | g), nu = 5), family = student(), data = d,
             importance = 200L)
  prop <- frmtmb:::imp_fit_proposal(fit)
  pl <- frmtmb:::imp_par_list(prop$template, fit$opt$par, prop$map)
  expect_equal(pl[["betad"]], fit$estimates[["betad"]])
  a <- as.matrix(prop$io$amat(pl))
  lme <- function(v) max(v) + log(mean(exp(v - max(v))))
  num <- -sum(apply(a, 1L, lme)) -
    sum(1 / 2 * log(2 * pi) + prop$plan$logdet_l)
  expect_equal(num, fit$opt$objective, tolerance = 1e-8)
  # the reported ESS describes the model: 0.50 when it read nu = sigma
  expect_gt(fit$importance$ess_min, 0.9)
})

test_that("a negative likelihood curvature is cut off, a positive one kept", {
  # the proposal adds a precision to this part, at any theta, so it has
  # to stay positive semi-definite; a log-concave likelihood must keep the
  # Laplace Hessian to the last bit
  pos <- matrix(c(2, 0.5, 0.5, 1), 2L)
  expect_identical(frmtmb:::imp_psd_part(pos), pos)
  neg <- matrix(c(1, 2, 2, 1), 2L)
  cut <- frmtmb:::imp_psd_part(neg)
  ev <- eigen(cut, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(ev[[1L]], 3)
  expect_lt(abs(ev[[2L]]) / ev[[1L]], 1e-12)
})

test_that("no count that works is not reported as one", {
  d <- imp_scalar_data()
  f <- bf(y ~ x + (1 | g))
  msg <- tryCatch(imp_with_displaced_rounds(
    1e9, frm(f, family = bernoulli(), data = d, importance = 8L)),
    error = conditionMessage)
  expect_match(msg, "at importance = 16, 32, 64, 128", fixed = TRUE)
  expect_match(msg, "No count above 128 was tried", fixed = TRUE)
  expect_false(grepl("settled at", msg, fixed = TRUE))
  expect_false(grepl("remedy", msg, fixed = TRUE))
  # at the top of its range it says so, by the rows it read
  local({
    testthat::local_mocked_bindings(imp_refuse_cells = 4000,
                                    .package = "frmtmb")
    msg <- tryCatch(imp_with_displaced_rounds(
      1e9, frm(f, family = bernoulli(), data = d, importance = 8L)),
      error = conditionMessage)
    expect_match(msg,
                 "on 400 rows it would tape more than 4e+03 rows times draws",
                 fixed = TRUE)
    expect_false(grepl("settled at", msg, fixed = TRUE))
  })
  # one doubling is all it read, and it claims no more than that
  local({
    testthat::local_mocked_bindings(imp_refuse_cells = 8000,
                                    .package = "frmtmb")
    msg <- tryCatch(imp_with_displaced_rounds(
      1e9, frm(f, family = bernoulli(), data = d, importance = 8L)),
      error = conditionMessage)
    expect_match(msg, "at importance = 16 either", fixed = TRUE)
    expect_false(grepl("remedy", msg, fixed = TRUE))
  })
})

test_that("a capped run has not settled, and a refused count is not named", {
  d <- imp_scalar_data()
  f <- bf(y ~ x + (1 | g))
  # above the wrapper's limit nothing rises, and a cap of one round keeps the
  # fit from settling: the old test read that as a count that works
  ctl <- frmtmb_control(importance_rounds = 1L)
  msg <- tryCatch(imp_with_displaced_rounds(
    64L, frm(f, family = bernoulli(), data = d, importance = 16L,
             control = ctl)),
    error = conditionMessage)
  expect_type(msg, "character")
  expect_false(grepl("settled at", msg, fixed = TRUE))
  expect_false(grepl("Use importance", msg, fixed = TRUE))
  expect_match(msg, "at importance = 32, 64, 128, 256 either", fixed = TRUE)
})

test_that("a slope near separation is never silent at few draws", {
  # 20 groups of 4, binary rows a slope nearly separates. Such a fit walks, and
  # where it lands is chaotic, so this asserts what must hold on either side: it
  # stops with the refusal, or it is capped and warns.
  set.seed(7L)
  g <- factor(rep(1:20, each = 4L))
  x <- rnorm(80L)
  eta <- x + rnorm(20L, 0, 2)[g] + rnorm(80L, 0, 0.1)
  d <- data.frame(y = as.integer(eta > 0),
                  x = x, g = g)
  ws <- character(0)
  r <- tryCatch(withCallingHandlers(
    frm(bf(y ~ x + (x | g)), family = bernoulli(), data = d, importance = 8L),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = conditionMessage)
  if (is.character(r)) {
    expect_match(r, "`importance = 8` did not settle", fixed = TRUE)
  } else {
    expect_true(r$importance$capped)
    expect_true(any(grepl("still moving", ws, fixed = TRUE)))
  }
})

test_that("an allocation failure in the ladder becomes the message", {
  real_plan <- frmtmb:::imp_plan
  real_opt <- frmtmb:::optimize_obj
  first <- TRUE
  d <- imp_scalar_data()
  msg <- tryCatch(testthat::with_mocked_bindings(
    frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d, importance = 16L),
    imp_plan = function(lap_obj, frame, lay, par, n_draw, seed) {
      if (n_draw > 32L) stop("std::bad_alloc", call. = FALSE)
      real_plan(lap_obj, frame, lay, par, n_draw, seed)
    },
    optimize_obj = function(obj, control, bounds = list(lower = -Inf,
                                                        upper = Inf),
                            par_units = NULL, verbose = 0L,
                            start_par = obj$par, tally = NULL) {
      o <- real_opt(obj, control, bounds, par_units, verbose, start_par, tally)
      if (!first) o$par[[1L]] <- o$par[[1L]] + 1
      first <<- FALSE
      o
    },
    .package = "frmtmb"), error = conditionMessage)
  expect_match(msg, "importance = 64 could not run: std::bad_alloc",
               fixed = TRUE)
  expect_match(msg, "at importance = 32 either", fixed = TRUE)
})

# ---------------------------------------------------------------------
# THE TAPE AGAINST AN EXACT VALUE, OFF THE ANCHOR
# ---------------------------------------------------------------------

test_that("the tape's curvature on a slope follows the moved mode", {
  # A step in a log-link intercept moves the mode by minus that step, so the
  # curvature at the moved mode stays put. Its partial derivative equals
  # itself, and a proposal that takes the partial halves this entry.
  set.seed(2102L)
  g <- factor(rep(1:10, each = 20L))
  x <- round(rnorm(200L), 3)
  a <- rnorm(10L)
  b <- 0.3 * a + sqrt(1 - 0.3^2) * rnorm(10L)
  eta <- 0.8 * a[g] + (0.3 + 0.4 * b[g]) * x + 0.5
  d <- data.frame(y = rpois(200L, exp(eta)),
                  x = x, g = g)
  p <- imp_parts(bf(y ~ x + (1 + x | g)) + poisson(), d)
  at <- p$opt$par
  h <- 1e-4
  e <- replace(numeric(length(at)), 1L, h)
  exact <- (p$lap$gr(at + e) - p$lap$gr(at - e))[[1L]] / (2 * h)
  for (s in 1:3) {
    plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, at, 200L, s)
    io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
    obj <- RTMB::MakeADFun(io$fn, frmtmb:::imp_template(p$tpl, "b"),
                           silent = TRUE)
    expect_lt(abs(obj$he(at)[1L, 1L] / exact - 1), 0.1)
  }
})

test_that("the objective two standard errors off the anchor is the integral", {
  # A wrong log-determinant is zero at the anchor and grows away from it, so
  # the anchor cannot see it: this reads an integrate() at displaced theta.
  d <- imp_textbook_data(124L)
  d$y <- rbinom(nrow(d), 10L, plogis(d$eta))
  p <- imp_parts(bf(y | trials(ntr) ~ x + (1 | g)) + binomial(), d)
  anchor <- p$opt$par
  lap <- frm(bf(y | trials(ntr) ~ x + (1 | g)), family = binomial(), data = d)
  se <- sqrt(diag(frmtmb:::sdr_of(lap)$cov.fixed))
  exact_nll <- function(par) {
    s1 <- exp(par[[3L]])
    -sum(vapply(split(seq_len(nrow(d)), d$g), function(r) {
      imp_log_integral(function(ui) {
        sum(dbinom(d$y[r], 10L, plogis(par[[1L]] + par[[2L]] * d$x[r] + ui),
                   log = TRUE)) + dnorm(ui, 0, s1, log = TRUE)
      }, s1)
    }, 0))
  }
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, anchor, 200L, 1L)
  io <- frmtmb:::build_importance_objective(p$frame, p$lay, p$gmap, plan)
  a0 <- as.matrix(io$amat(frmtmb:::imp_par_list(p$tpl, anchor)))
  io$set_offset(apply(a0, 1L, max))
  obj <- RTMB::MakeADFun(io$fn, frmtmb:::imp_template(p$tpl, "b"),
                         silent = TRUE)
  for (sgn in list(c(2, 2, 2), c(-2, 2, -2), c(0, 0, 2), c(0, 0, -2))) {
    at <- anchor + sgn * se
    mcse <- frmtmb:::imp_ess(
      as.matrix(io$amat(frmtmb:::imp_par_list(p$tpl, at))),
      plan[["n_draw"]])$mcse
    expect_lt(abs(obj$fn(at) - exact_nll(at)), 4 * mcse)
  }
})

test_that("a fit does not read a column's center or scale", {
  # the same model on an affine map of its column: the slope's SE scales with
  # it, and the log-likelihood and theta's SE do not move. A step in the
  # coefficient's own units gave SE ratios of 1.47, capped rounds and refusals
  # on a calendar year.
  set.seed(7L)
  g <- factor(rep(1:15, each = 10L))
  z <- rnorm(150L)
  y <- rpois(150L, exp(0.5 + 0.3 * z + rnorm(15L, 0, 0.8)[g]))
  fit_on <- function(a, b) {
    d <- data.frame(y = y, w = a * z + b, g = g)
    allow_warnings(
      f <- frm(bf(y ~ w + (1 | g)), family = poisson(), data = d,
               importance = 200L),
      "false convergence")
    f
  }
  ref <- fit_on(1, 0)
  se_ref <- sqrt(diag(frmtmb:::sdr_of(ref)$cov.fixed))
  for (m in list(c(1, 2000), c(1000, 0), c(0.001, 0), c(100, 1e4))) {
    f <- fit_on(m[[1L]], m[[2L]])
    expect_false(f$importance$capped)
    expect_lt(abs(as.numeric(logLik(f)) - as.numeric(logLik(ref))) /
                ref$importance$mcse, 0.05)
    se <- sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed))
    expect_lt(abs(se[[2L]] * m[[1L]] / se_ref[[2L]] - 1), 0.01)
    expect_lt(abs(se[[3L]] / se_ref[[3L]] - 1), 0.01)
    expect_lt(abs(f$opt$par[[2L]] * m[[1L]] - ref$opt$par[[2L]]) / se_ref[[2L]],
              0.01)
  }
})

test_that("a dispersion's column is read the same way", {
  set.seed(9L)
  g <- factor(rep(1:20, each = 8L))
  z <- rnorm(160L)
  w <- rnorm(160L)
  y <- rnbinom(160L, mu = exp(1 + 0.3 * z + rnorm(20L, 0, 0.6)[g]),
               size = exp(0.7 + 0.3 * w))
  fit_on <- function(b) {
    d <- data.frame(y = y, z = z, w = w + b, g = g)
    frm(bf(y ~ z + (1 | g), shape ~ w), family = negbinomial(), data = d,
        importance = 200L)
  }
  ref <- fit_on(0)
  f <- fit_on(2000)
  expect_false(f$importance$capped)
  expect_lt(abs(as.numeric(logLik(f)) - as.numeric(logLik(ref))) /
              ref$importance$mcse, 0.05)
  se_ref <- sqrt(diag(frmtmb:::sdr_of(ref)$cov.fixed))
  se <- sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed))
  j <- which(names(ref$opt$par) == "betad")[[2L]]
  expect_lt(abs(se[[j]] / se_ref[[j]] - 1), 0.02)
  # the fit read the column, not luck: its step is the column's
  steps <- frmtmb:::imp_outer_steps(f$frame, f$opt$par)
  expect_lt(steps[[j]] / frmtmb:::imp_cross_step, 1e-3)
})

test_that("a negative curvature of the likelihood is cut, on a real fit", {
  # Student rows that disagree make the likelihood's curvature at the mode
  # negative: the mode sits between them, where each density is convex. The
  # proposal then rests on the prior's precision alone. abs() in place of the
  # cut gives a narrower Gaussian at a group that needs a wider one.
  set.seed(2L)
  g <- factor(rep(1:15, each = 2L))
  x <- rnorm(30L)
  y <- 0.3 + 0.5 * x + rnorm(15L, 0, 0.7)[g] + 0.3 * rt(30L, 3L)
  y[c(1L, 7L, 13L)] <- y[c(1L, 7L, 13L)] + 12
  d <- data.frame(y = y, x = x, g = g)
  p <- imp_parts(bf(y ~ x + (1 | g), nu = 2) + student(), d)
  at <- p$opt$par
  plan <- frmtmb:::imp_plan(p$lap, p$frame, p$lay, at, 100L, 1L)
  cut <- which(plan$hlik[1L, 1L, ] == 0)
  expect_gt(length(cut), 0L)
  # at a cut group the proposal's variance is the prior's
  k <- cut[[1L]]
  ii <- plan$idx[1L, k]
  sdv <- exp(plan$logdet_l[[k]])
  expect_equal(sdv^2 * plan$prec[1L, 1L], 1)
  expect_equal(unname((plan$u[ii, 1L] - plan$uhat[ii]) / sdv),
               unname(plan$z[ii, 1L]))
})

test_that("a prior enters the corrected objective", {
  # It used to stop with \"Report this model\": the pin compared the joint
  # density, with the prior, against pieces without it.
  d <- imp_textbook_data(124L)
  d$y <- rbinom(nrow(d), 10L, plogis(d$eta))
  f <- bf(y | trials(ntr) ~ x + (1 | g))
  for (pr in list(set_prior("normal(0, 1)", class = "b"),
                  set_prior("exponential(1)", class = "sd"))) {
    qd <- frm(f, family = binomial(), data = d, quadrature = TRUE, prior = pr)
    fit <- frm(f, family = binomial(), data = d, importance = 200L, prior = pr)
    expect_lt(abs(fit$opt$objective - qd$opt$objective),
              3 * fit$importance$mcse)
    se <- sqrt(diag(frmtmb:::sdr_of(qd)$cov.fixed))
    expect_lt(max(abs(fit$opt$par - qd$opt$par) / se), 0.05)
    # and it is a prior: without it the fit is another fit
    plain <- frm(f, family = binomial(), data = d, importance = 200L)
    expect_gt(max(abs(fit$opt$par - plain$opt$par) / se), 1e-3)
  }
  expect_identical(frm_compat("importance", "prior")$status, "works")
})

test_that("the rounds reach the exact optimum where Laplace is far from it", {
  # Laplace puts the intercept at -8.9 and the exact optimum is -2.5, so one
  # proposal anchored there is not enough: the first round lands 1.3 standard
  # errors away, and the rounds after it come to 0.1.
  set.seed(1003L)
  g <- factor(rep(seq_len(40L), each = 3L))
  x <- round(rnorm(120L), 3)
  u <- rnorm(40L, 0, 4)
  d <- data.frame(y = rbinom(120L, 1, plogis(-0.5 + 0.5 * x + u[g])), x = x,
                  g = g)
  f <- bf(y ~ x + (1 | g))
  qd <- frm(f, family = bernoulli(), data = d, quadrature = TRUE)
  se <- sqrt(diag(frmtmb:::sdr_of(qd)$cov.fixed))
  allow_warnings(
    fit <- frm(f, family = bernoulli(), data = d, importance = 200L),
    c("covers", "still moving"))
  expect_lt(max(abs(fit$opt$par - qd$opt$par) / se), 0.5)
})

test_that("a count whose worst group is below the floor is not named", {
  # On a slope near separation, 128 draws passed at a worst group of 0.012 with
  # an error of 1.87 nats, and 256 to 1024 capped. A floor no group clears
  # stands in for that here: the same wrapper settles, and names nothing.
  d <- imp_scalar_data()
  f <- bf(y ~ x + (1 | g))
  msg <- tryCatch(imp_with_displaced_rounds(
    64L, frm(f, family = bernoulli(), data = d, importance = 16L,
             control = frmtmb_control(importance_ess = 0.9999))),
    error = conditionMessage)
  expect_type(msg, "character")
  expect_false(grepl("Use importance", msg, fixed = TRUE))
  expect_false(grepl("settled at", msg, fixed = TRUE))
  expect_match(msg, "above an effective sample size of 0.9999", fixed = TRUE)
})

# Ordinal rows on a column whose center the first threshold takes.
imp_ord_data <- function(seed, m) {
  set.seed(seed)
  g <- factor(rep(seq_len(15L), each = 10L))
  z <- rnorm(150L)
  lat <- 0.7 * z + rnorm(15L, 0, 0.7)[g] + rlogis(150L)
  data.frame(y = ordered(as.integer(cut(lat, c(-Inf, -1, 0.5, Inf)))),
             w = 5 * z + m, g = g)
}

test_that("a threshold does not step by its own size", {
  # at m = 1e4 the first threshold is about 1600. A step of its own size moved
  # the predictor by 0.16 and read a slope SE of 4.56 times the centered fit,
  # with no warning; at 3e4 the SE was NaN and the rounds capped.
  f <- bf(y ~ w + (1 | g))
  ref <- frm(f, family = cumulative(), data = imp_ord_data(5L, 0),
             importance = 300L)
  se_ref <- sqrt(diag(frmtmb:::sdr_of(ref)$cov.fixed))
  j <- which(names(ref$opt$par) == "beta")
  for (m in c(1e4, 3e4)) {
    allow_warnings(
      fit <- frm(f, family = cumulative(), data = imp_ord_data(5L, m),
                 importance = 300L),
      "false convergence")
    expect_false(fit$importance$capped)
    expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))) /
                ref$importance$mcse, 0.05)
    se <- sqrt(diag(frmtmb:::sdr_of(fit)$cov.fixed))
    expect_lt(abs(se[[j]] / se_ref[[j]] - 1), 0.01)
    expect_lt(abs(fit$opt$par[[j]] - ref$opt$par[[j]]) / se_ref[[j]], 0.01)
    # and no step reads a parameter's own size
    steps <- frmtmb:::imp_outer_steps(fit$frame, fit$opt$par)
    tau <- names(fit$opt$par) == "tau_raw"
    expect_true(any(abs(fit$opt$par[tau]) > 100))
    expect_true(all(steps[tau] == frmtmb:::imp_cross_step))
  }
})

test_that("a sd far from one reads its own precision", {
  # A precision read as a difference at a unit vector rounds once 1 / sd^2 nears
  # the normalizer's eps: 0.29 wrong at sd 1e7, and 0 at 1e8. The pin then
  # stopped
  # with a request to report the model.
  set.seed(3L)
  g <- factor(rep(seq_len(15L), each = 10L))
  z <- rnorm(150L)
  u <- rnorm(15L, 0, 0.7)[g]
  y0 <- 0.3 + 0.5 * z + u + rnorm(150L, 0, 0.5)
  yt <- 0.3 + 0.5 * z + u + 0.5 * rt(150L, 5L)
  ref <- frm(bf(y ~ z + (1 | g), nu = 4), family = student(),
             data = data.frame(y = yt, z = z, g = g), importance = 200L)
  for (k in c(1e6, 1e7)) {
    d <- data.frame(y = y0 * k, z = z, g = g)
    lap <- frm(bf(y ~ z + (1 | g)), family = gaussian(), data = d)
    fit <- frm(bf(y ~ z + (1 | g)), family = gaussian(), data = d,
               importance = 64L)
    expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(lap))) /
                abs(as.numeric(logLik(lap))), 1e-8)
    d$y <- yt * k
    st <- frm(bf(y ~ z + (1 | g), nu = 4), family = student(), data = d,
              importance = 200L)
    expect_lt(abs(as.numeric(logLik(st)) + 150 * log(k) -
                    as.numeric(logLik(ref))),
              3 * ref$importance$mcse)
  }
  # the reading itself, against the sd it was given
  lay <- list(q = 1L, offset = 0L, blocks = lap$frame[["re_blocks"]])
  for (th in log(c(1e-3, 1, 1e8))) {
    tt <- frmtmb:::imp_prior_scales(lay, th)
    tm <- frmtmb:::imp_prior_terms(lay, tt)(th)
    expect_lt(abs(-2 * tm$diag[[1L]] * exp(2 * th) - 1), 1e-8)
  }
})

test_that("the rounds start where the Laplace route ends", {
  # A column's center makes a false collapse. The Laplace route leaves it, and a
  # correction anchored there has every weight equal and nothing left to do: at
  # data seed 8 the fit stayed at log sd -13.5, 0.79 nats below the centered
  # fit.
  f <- bf(y ~ w + (1 | g))
  d <- imp_ord_data(8L, 2000)
  ref <- frm(f, family = cumulative(), data = imp_ord_data(8L, 0),
             importance = 200L)
  allow_warnings(
    fit <- frm(f, family = cumulative(), data = d, importance = 200L),
    "false convergence")
  th <- which(names(fit$opt$par) == "theta")
  se <- sqrt(diag(frmtmb:::sdr_of(ref)$cov.fixed))
  expect_lt(abs(fit$opt$par[[th]] - ref$opt$par[[th]]) / se[[th]], 0.05)
  expect_lt(abs(as.numeric(logLik(fit)) - as.numeric(logLik(ref))) /
              ref$importance$mcse, 0.5)
})

test_that("a predictor in small response units steps in those units", {
  # a step of 1e-4 is 100 of this response's scale. Without the unit, this fit
  # stopped on 2 of 4 data seeds and read SE ratios to 1.007 on the others.
  for (s in 3:4) {
    set.seed(s)
    g <- factor(rep(seq_len(15L), each = 10L))
    z <- rnorm(150L)
    y0 <- 0.3 + 0.5 * z + rnorm(15L, 0, 0.7)[g] + 0.5 * rt(150L, 5L)
    f <- bf(y ~ z + (1 | g), nu = 4)
    ref <- frm(f, family = student(), data = data.frame(y = y0, z = z, g = g),
               importance = 200L)
    fit <- frm(f, family = student(),
               data = data.frame(y = y0 * 1e-6, z = z, g = g),
               importance = 200L)
    expect_lt(abs(as.numeric(logLik(fit)) + 150 * log(1e-6) -
                    as.numeric(logLik(ref))), 3 * ref$importance$mcse)
    se <- sqrt(diag(frmtmb:::sdr_of(fit)$cov.fixed))
    se_ref <- sqrt(diag(frmtmb:::sdr_of(ref)$cov.fixed))
    expect_lt(abs(se[[2L]] / 1e-6 / se_ref[[2L]] - 1), 0.001)
  }
})
