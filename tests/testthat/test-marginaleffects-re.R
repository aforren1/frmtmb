# marginaleffects on fits with group-level terms, smooths and fields:
# get_coef() carries the group-level effects (and the theta a field's
# coefficients read), and get_vcov() their joint covariance with the
# coefficients, so the delta method holds the uncertainty fitted() holds.
# Through 0.79.1 get_coef() held the coefficients alone: on 168 mixed fits
# avg_predictions() was off by more than 25 percent on 77 percent of them
# (dev/review-predread-r2-mfxrate.R; dev/mfxse-findings.md has the after).

mre_data <- function(seed = 7411, J = 12L, m = 20L) {
  set.seed(seed)
  n <- J * m
  d <- data.frame(x = stats::rnorm(n), f = factor(sample(c("a", "b"), n,
                                                         TRUE)),
                  g = factor(rep(seq_len(J), each = m)))
  u0 <- stats::rnorm(J, 0, 0.8)
  u1 <- stats::rnorm(J, 0, 0.5)
  d$cnt <- stats::rpois(n, exp(0.5 + (0.4 + u1[d$g]) * d$x +
                                 0.3 * (d$f == "b") + u0[d$g]))
  d$y <- 1 + 0.5 * d$x + u0[d$g] + stats::rnorm(n, 0, 0.6)
  d
}

# the coefficient block alone: the standard error as 0.79.1 formed it, with
# the group-level effects held at their modes
mre_coef_only <- function(fit) {
  nm <- interop_coef_names(fit)
  marginaleffects::get_vcov(fit)[nm, nm]
}
mre_held <- "covers the coefficients of get_coef() alone"

# The delta method by hand: the Jacobian of a mean over rows from
# frm_lp_basis(), which knows nothing of get_coef() or get_vcov(). `inv` maps
# eta to d mu / d eta; with `nd_lo`, the mean is of a difference over `h`.
mre_avg_se <- function(fit, nd_hi, nd_lo = NULL, h = 1, inv = exp) {
  bh <- frm_lp_basis(fit, newdata = nd_hi)
  Jh <- as.matrix(bh$A) * inv(bh$eta)
  J <- colMeans(Jh)
  if (!is.null(nd_lo)) {
    bl <- frm_lp_basis(fit, newdata = nd_lo)
    expect_identical(bl$coef_pos, bh$coef_pos)
    J <- (J - colMeans(as.matrix(bl$A) * inv(bl$eta))) / h
  }
  sqrt(sum(J * (bh$V %*% J)))
}

test_that("predictions() on a mixed fit reads fitted()'s Est.Error", {
  skip_if_not_installed("marginaleffects")
  d <- mre_data()
  fit <- frm(cnt ~ x + f + (1 + x | g), family = poisson(), data = d)
  nd <- d[c(1, 50, 120, 200), ]
  p <- marginaleffects::predictions(fit, newdata = nd)
  ft <- fitted(fit, newdata = nd)[, "Est.Error"]
  old <- allow_warnings(
    marginaleffects::predictions(fit, newdata = nd,
                                 vcov = mre_coef_only(fit)),
    mre_held, require = mre_held)
  # the yardstick is how far the 0.79.1 convention is from fitted(): the
  # agreement must be three orders of magnitude closer than its worst row
  gap_old <- max(abs(old$std.error / ft - 1))
  expect_gt(gap_old, 0.2)
  expect_lt(max(abs(p$std.error / ft - 1)), 1e-3 * gap_old)
  # and the estimate is fitted()'s, as it was
  expect_equal(p$estimate, unname(fitted(fit, newdata = nd)[, "Estimate"]))
})

test_that("an average and a slope hold the group-level uncertainty", {
  skip_if_not_installed("marginaleffects")
  d <- mre_data()
  fit <- frm(cnt ~ x + f + (1 + x | g), family = poisson(), data = d)
  ap <- marginaleffects::avg_predictions(fit)
  ref <- mre_avg_se(fit, d)
  h <- 1e-4
  dh <- transform(d, x = x + h / 2)
  dl <- transform(d, x = x - h / 2)
  sl <- marginaleffects::avg_slopes(fit, variables = "x")
  ref_sl <- mre_avg_se(fit, dh, dl, h)
  old_ap <- allow_warnings(marginaleffects::avg_predictions(
    fit, vcov = mre_coef_only(fit)), mre_held)
  old_sl <- allow_warnings(marginaleffects::avg_slopes(
    fit, variables = "x", vcov = mre_coef_only(fit)), mre_held)
  expect_lt(abs(ap$std.error - ref) / abs(old_ap$std.error - ref), 1e-2)
  expect_lt(abs(sl$std.error - ref_sl) / abs(old_sl$std.error - ref_sl),
            1e-2)
  # 0.79.1 was too wide on both, by more than a quarter on the average
  expect_gt(old_ap$std.error / ref, 1.25)
})

test_that("re_formula = NA reads the coefficients alone, as fitted() does", {
  skip_if_not_installed("marginaleffects")
  d <- mre_data()
  fit <- frm(cnt ~ x + f + (1 + x | g), family = poisson(), data = d)
  nd <- d[c(1, 50, 120), ]
  w <- "not known to be supported"
  cf <- marginaleffects::get_coef(fit, re_formula = NA)
  expect_false(any(startsWith(names(cf), "b.")))
  expect_true(any(startsWith(names(marginaleffects::get_coef(fit)), "b.")))
  p <- allow_warnings(marginaleffects::predictions(fit, newdata = nd,
                                                   re_formula = NA), w)
  ft <- fitted(fit, newdata = nd, re_formula = NA)[, "Est.Error"]
  expect_equal(p$std.error, unname(ft), tolerance = 1e-6)
})

test_that("hypotheses() on the coefficients keeps their rows and errors", {
  skip_if_not_installed("marginaleffects")
  d <- mre_data()
  fit <- frm(cnt ~ x + f + (1 + x | g), family = poisson(), data = d)
  h <- marginaleffects::hypotheses(fit)
  fe <- fixef(fit)
  expect_identical(nrow(h), nrow(fe))
  expect_equal(as.numeric(h$std.error), unname(fe[, "Est.Error"]),
               tolerance = 1e-6)
  w <- "positional indices"
  h2 <- allow_warnings(marginaleffects::hypotheses(fit, "b2 = 0"), w)
  expect_equal(as.numeric(h2$estimate), unname(fe["x", "Estimate"]))
  expect_equal(as.numeric(h2$std.error), unname(fe["x", "Est.Error"]),
               tolerance = 1e-6)
  # the default joint test is of the coefficients: marginaleffects' joint
  # test calls get_coef() with no argument, and with the group-level effects
  # in it the test was of 24 modes as well (F 328 on 27 df against 15.4 on
  # 3; dev/reviews/2026-11-06-mfxse.md, M1). The reference is the Wald
  # statistic from vcov() by hand
  hj <- marginaleffects::hypotheses(fit, joint = TRUE)
  b <- fe[, "Estimate"]
  V <- vcov(fit)
  expect_identical(as.integer(hj$df1), length(b))
  expect_equal(as.numeric(hj$statistic),
               drop(t(b) %*% solve(V, b)) / length(b), tolerance = 1e-6)
  # a pattern that is not a coefficient name matches no group-level name
  expect_error(marginaleffects::hypotheses(fit, joint = "^b[.]"),
               "does not match any estimate name")
})

# A function as `hypothesis` gets the model and is differentiated as any
# quantity is, so the effects stay in: holding them there read 1.774 where
# avg_predictions() alone reads 0.143 (dev/reviews/2026-11-06-mfxse.md,
# short check, S1)
test_that("hypotheses(fit, <function>) holds the group-level effects", {
  skip_if_not_installed("marginaleffects")
  d <- mre_data()
  fit <- frm(cnt ~ x + f + (1 + x | g), family = poisson(), data = d)
  ap <- marginaleffects::avg_predictions(fit)
  hf <- marginaleffects::hypotheses(fit, hypothesis = function(x) {
    data.frame(term = "avg",
               estimate = marginaleffects::avg_predictions(x)$estimate)
  })
  old <- allow_warnings(marginaleffects::avg_predictions(
    fit, vcov = mre_coef_only(fit)), mre_held)
  gap_old <- abs(old$std.error / ap$std.error - 1)
  expect_gt(gap_old, 1)
  expect_lt(abs(as.numeric(hf$std.error) / ap$std.error - 1),
            1e-3 * gap_old)
})

test_that("set_coef() moves the rows of the group it sets, and no other", {
  skip_if_not_installed("marginaleffects")
  d <- mre_data()
  fit <- frm(cnt ~ x + f + (1 + x | g), family = poisson(), data = d)
  cf <- marginaleffects::get_coef(fit)
  k <- which(names(cf) == "b.1 + x | g.3.(Intercept)")
  expect_length(k, 1L)
  cf2 <- cf
  cf2[k] <- cf2[k] + 0.5
  f2 <- marginaleffects::set_coef(fit, cf2)
  expect_equal(unname(marginaleffects::get_coef(f2)), unname(cf2))
  p0 <- frm_linpred(fit, newdata = d)
  p2 <- frm_linpred(f2, newdata = d)
  moved <- abs(p2 - p0) > 0
  expect_true(all(moved[d$g == "3"]))
  expect_false(any(moved[d$g != "3"]))
})

test_that("a smooth and a gp() field are in get_coef(), with no warning", {
  skip_if_not_installed("marginaleffects")
  set.seed(7411)
  n <- 200
  d <- data.frame(x = stats::rnorm(n))
  d$y <- sin(2 * d$x) + stats::rnorm(n, 0, 0.4)
  nd <- data.frame(x = c(-1, 0, 1, 2.5))
  for (f in list(frm(y ~ s(x), data = d), frm(y ~ gp(x, k = 15), data = d))) {
    p <- marginaleffects::predictions(f, newdata = nd)
    ft <- fitted(f, newdata = nd)[, "Est.Error"]
    # 0.79.1 read 6.9, 0.51, 6.5 and 6.8 times fitted()'s on s(x)
    expect_equal(p$std.error, unname(ft), tolerance = 1e-5)
  }
})

test_that("avg_slopes() on an exact gp() field holds the field", {
  skip_if_not_installed("marginaleffects")
  set.seed(8150)
  m <- 90
  d <- data.frame(x = stats::runif(m, 0, 6), z = stats::runif(m, 0, 6))
  d$y <- sin(d$x) + 0.3 * d$z + stats::rnorm(m, sd = 0.3)
  fit <- frm(y ~ z + gp(x), data = d)
  sl <- marginaleffects::avg_slopes(fit, variables = "x")
  h <- 1e-4
  ref <- mre_avg_se(fit, transform(d, x = x + h / 2),
                    transform(d, x = x - h / 2), h, inv = function(e) 1)
  # 0.79.1 read 1e-6 or so here: get_coef() held no parameter of the field
  expect_gt(sl$std.error, 0.5 * ref)
  expect_lt(abs(sl$std.error / ref - 1), 1e-2)
})

test_that("what no parameter holds is still said, by name", {
  skip_if_not_installed("marginaleffects")
  d <- mre_data()
  fit <- frm(cnt ~ x + f + (1 | g), family = poisson(), data = d)
  nd <- data.frame(x = c(0, 1), f = "a", g = factor(c("3", "new")))
  # a level the fit did not see: fitted() adds the block's variance there
  w <- "1 read a grouping level the fit did not see"
  p <- allow_warnings(marginaleffects::predictions(
    fit, newdata = nd, allow_new_levels = TRUE),
    c(w, "not known to be supported"), require = w)
  ft <- fitted(fit, newdata = nd, allow_new_levels = TRUE)[, "Est.Error"]
  expect_equal(p$std.error[1L], unname(ft[1L]), tolerance = 1e-5)
  expect_lt(p$std.error[2L], ft[2L])
  # a slope and a contrast say it too. comparisons() makes its estimate pass
  # after get_vcov(), inside a handler that muffles every warning, and the
  # check that warned there said nothing (dev/reviews/2026-11-06-mfxse.md,
  # M2); it now warns from the first prediction of the Jacobian
  nd3 <- data.frame(x = c(0, 1, 2), f = factor("a", levels = c("a", "b")),
                    g = factor(c("3", "new", "new2")))
  w2 <- "2 read a grouping level the fit did not see"
  allow_warnings(marginaleffects::avg_slopes(
    fit, newdata = nd3, variables = "x", allow_new_levels = TRUE),
    c(w2, "not known to be supported"), require = w2)
  allow_warnings(marginaleffects::avg_comparisons(
    fit, newdata = nd3, variables = "f", allow_new_levels = TRUE),
    c(w2, "not known to be supported"), require = w2)
  # a user's covariance holds no new level's variance either, and says so
  for (v in list(~ g, vcov(fit, full = TRUE))) {
    allow_warnings(marginaleffects::predictions(
      fit, newdata = nd3, vcov = v, allow_new_levels = TRUE),
      c(w2, "not known to be supported"), require = w2)
  }
  # the absent case: seen levels only, no warning of any kind
  expect_no_warning(p0 <- marginaleffects::predictions(fit,
                                                       newdata = d[1:3, ]))
  expect_true(all(is.finite(p0$std.error)))
  # a matrix over the coefficients alone holds the effects at their modes
  V <- mre_coef_only(fit)
  allow_warnings(marginaleffects::predictions(fit, newdata = d[1:3, ],
                                              vcov = V),
                 mre_held, require = mre_held)
  # a fit that integrates the effects out has no covariance for them
  fq <- frm(cnt ~ x + f + (1 | g), family = poisson(), data = d,
            quadrature = TRUE)
  wq <- "integrates out"
  allow_warnings(marginaleffects::predictions(fq, newdata = d[1:3, ]),
                 c(wq, "conditional modes"), require = wq)
})

test_that("an exact gp() far from its data warns; a slope near it does not", {
  skip_if_not_installed("marginaleffects")
  set.seed(8150)
  m <- 90
  d <- data.frame(x = stats::runif(m, 0, 6), z = stats::runif(m, 0, 6))
  d$y <- sin(d$x) + 0.3 * d$z + stats::rnorm(m, sd = 0.3)
  fit <- frm(y ~ z + gp(x), data = d)
  w <- "1 read an exact gp() at a position the fit did not see"
  nd <- data.frame(x = c(d$x[5], 3.33, 12), z = 3)
  p <- allow_warnings(marginaleffects::predictions(fit, newdata = nd), w,
                      require = w)
  ft <- fitted(fit, newdata = nd)[, "Est.Error"]
  expect_lt(p$std.error[3L], 0.8 * ft[3L])
  # the rows at and between the data agree; the slope reads x + eps, an
  # unseen position with a kriging variance of order eps^2, and says nothing
  expect_equal(p$std.error[1:2], unname(ft[1:2]), tolerance = 1e-4)
  sl <- marginaleffects::avg_slopes(fit, variables = "x")
  expect_true(is.finite(sl$std.error))
})

test_that("a new level and an unseen gp() position are each counted", {
  skip_if_not_installed("marginaleffects")
  set.seed(8150)
  m <- 90
  d <- data.frame(x = stats::runif(m, 0, 6),
                  g = factor(sample(1:6, m, TRUE)))
  d$y <- sin(d$x) + stats::rnorm(6, 0, 0.5)[d$g] + stats::rnorm(m, sd = 0.3)
  fit <- allow_warnings(frm(y ~ gp(x) + (1 | g), data = d), "boundary")
  nd <- data.frame(x = c(d$x[1], 3.33, 9), g = factor(c("1", "new", "2")))
  w <- c("1 read a grouping level the fit did not see",
         "1 read an exact gp() at a position the fit did not see")
  allow_warnings(marginaleffects::predictions(fit, newdata = nd,
                                              allow_new_levels = TRUE),
                 c(w, "not known to be supported"), require = w)
})

# A fit with two blocks: every block's effects are in get_coef(). A build
# that kept the first block alone read avg_predictions() 2.73 times the
# reference on such a fit and passed every other test
# (dev/reviews/2026-11-06-mfxse.md, M3, mutant oneblock)
test_that("crossed blocks, and a smooth beside a group, are all held", {
  skip_if_not_installed("marginaleffects")
  set.seed(9102)
  n <- 240
  d <- data.frame(x = stats::rnorm(n), g = factor(sample(1:12, n, TRUE)),
                  h = factor(sample(1:8, n, TRUE)))
  d$cnt <- stats::rpois(n, exp(0.4 + 0.3 * d$x +
                                 stats::rnorm(12, 0, 0.6)[d$g] +
                                 stats::rnorm(8, 0, 0.6)[d$h]))
  d$y <- sin(2 * d$x) + stats::rnorm(12, 0, 0.6)[d$g] +
    stats::rnorm(n, 0, 0.4)
  fits <- list(frm(cnt ~ x + (1 | g) + (1 | h), family = poisson(), data = d),
               frm(y ~ s(x) + (1 | g), data = d))
  invs <- list(exp, function(e) 1)
  nd <- d[c(2, 60, 130, 200), ]
  for (k in seq_along(fits)) {
    fit <- fits[[k]]
    p <- marginaleffects::predictions(fit, newdata = nd)
    ft <- fitted(fit, newdata = nd)[, "Est.Error"]
    old <- allow_warnings(marginaleffects::predictions(
      fit, newdata = nd, vcov = mre_coef_only(fit)), mre_held)
    gap_old <- max(abs(old$std.error / ft - 1))
    expect_gt(gap_old, 0.05)
    expect_lt(max(abs(p$std.error / ft - 1)), 1e-3 * gap_old)
    ap <- marginaleffects::avg_predictions(fit)
    ref <- mre_avg_se(fit, d, inv = invs[[k]])
    old_ap <- allow_warnings(marginaleffects::avg_predictions(
      fit, vcov = mre_coef_only(fit)), mre_held)
    expect_lt(abs(ap$std.error / ref - 1),
              1e-2 * abs(old_ap$std.error / ref - 1))
  }
})

test_that("vcov = carries an outer covariance into the group-level effects", {
  skip_if_not_installed("marginaleffects")
  d <- mre_data()
  fit <- frm(y ~ x + (1 | g), data = d)
  nd <- d[c(1, 50, 120), ]
  p <- marginaleffects::predictions(fit, newdata = nd)
  # the model's own outer covariance, carried, is the default
  pf <- marginaleffects::predictions(fit, newdata = nd,
                                     vcov = vcov(fit, full = TRUE))
  expect_equal(pf$std.error, p$std.error, tolerance = 1e-6)
  # a multiple of get_vcov() over all of its rows doubles every variance
  p2 <- marginaleffects::predictions(
    fit, newdata = nd, vcov = function(x) 2 * marginaleffects::get_vcov(x))
  expect_equal(p2$std.error / p$std.error, rep(sqrt(2), 3), tolerance = 1e-6)
  # a cluster sandwich over the groups is carried too, and differs
  ps <- marginaleffects::predictions(fit, newdata = nd, vcov = ~ g)
  expect_true(all(is.finite(ps$std.error)))
  expect_false(isTRUE(all.equal(ps$std.error, p$std.error)))
  # the carry forms a dense block over the inner rows, so it has the dense
  # limit of a latent block, and says so when it holds the effects instead
  withr::local_options(frmtmb.latres_dense_rows = 5L)
  w <- "could not be carried into the group-level effects"
  allow_warnings(marginaleffects::predictions(fit, newdata = nd, vcov = ~ g),
                 w, require = w)
})

# The sandwich reaches the effects' block, not the coefficients' alone. A
# build that took the effects' block from the model passed the test above
# (dev/reviews/2026-11-06-mfxse.md, M3, mutant sandown). The reference is
# the joint covariance by hand, Qrr^-1 + J Vf J' with Vf the sandwich, read
# through frm_lp_basis()'s design, on data whose noise grows with x, where
# the sandwich moves the effects' block.
test_that("a cluster sandwich is carried into the effects' block", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("sandwich")
  set.seed(9301)
  J <- 30L
  m <- 10L
  d <- data.frame(x = stats::rnorm(J * m),
                  g = factor(rep(seq_len(J), each = m)))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(J, 0, 0.7)[d$g] +
    stats::rnorm(J * m, 0, exp(0.6 * d$x) * 0.5)
  fit <- frm(y ~ x + (1 | g), data = d)
  nd <- d[c(1, 101, 201), ]
  p0 <- marginaleffects::predictions(fit, newdata = nd)
  ps <- marginaleffects::predictions(fit, newdata = nd, vcov = ~ g)
  Q <- joint_precision(fit)
  r <- fit$obj$env$random
  f <- setdiff(seq_len(nrow(Q)), r)
  Vf <- vcov_cluster(fit, ~ g, type = "CR1", full = TRUE)
  W <- solve(as.matrix(Q[r, r]))
  Jb <- -W %*% as.matrix(Q[r, f])
  V <- matrix(0, nrow(Q), nrow(Q))
  V[r, r] <- W + Jb %*% Vf %*% t(Jb)
  V[r, f] <- Jb %*% Vf
  V[f, r] <- t(V[r, f])
  V[f, f] <- Vf
  lb <- frm_lp_basis(fit, newdata = nd)
  A <- as.matrix(lb$A)
  Vc <- V[lb$coef_pos, lb$coef_pos]
  ref <- sqrt(rowSums((A %*% Vc) * A))
  moved <- max(abs(ref / p0$std.error - 1))
  expect_gt(moved, 1e-3)
  expect_lt(max(abs(ps$std.error / ref - 1)), 0.05 * moved)
})
