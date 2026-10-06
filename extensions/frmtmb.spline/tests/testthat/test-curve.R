## The curve covariance, and the check that licenses it.
##
## frmtmb exports no route to the joint covariance of a grid prediction,
## so this package rebuilds it. Every assertion here is about that
## rebuild being the same object core uses, not about it being plausible.

sp_curve_fit <- function(n = 250, seed = 7, k = 10) {
  set.seed(seed)
  d <- data.frame(x = sort(stats::runif(n)))
  d$y <- 2 * sin(pi * d$x) + 0.6 * d$x + stats::rnorm(n, 0, 0.4)
  list(d = d,
       fit = frmtmb::frm(frmtmb::bf(y ~ s(x, k = k)),
                         family = stats::gaussian(), data = d))
}

test_that("the assembled covariance reproduces frm_linpred(se.fit) exactly", {
  o <- sp_curve_fit()
  g <- data.frame(x = seq(0, 1, length.out = 40))
  cv <- frm_curve(o$fit, newdata = g, nsim = 2000, seed = 1)
  ck <- attr(cv, "check")
  # machine precision, not a tolerance chosen to pass
  expect_lt(ck$cov_rel_error, 1e-10)
  p <- frm_linpred(o$fit, newdata = g, type = "link", se.fit = TRUE,
                      re_formula = NA)
  expect_equal(cv$.estimate, as.numeric(p$fit), tolerance = 1e-12)
  expect_equal(cv$.se, as.numeric(p$se.fit), tolerance = 1e-10)
  # and the diagonal of the returned matrix is that same standard error
  expect_equal(sqrt(diag(attr(cv, "Sigma"))), cv$.se, tolerance = 1e-12)
})

test_that("the design costs one frm_linpred() call whatever its width", {
  # Up to frmtmb 0.51.0 this package rebuilt the design by unit
  # perturbation and s(x, k = 10) cost 11 frm_linpred() calls: 10 live
  # columns plus sigma's intercept, which cannot move the mu curve and
  # whose column came back all zero. frm_lp_basis() returns the design
  # core already had, so the only frm_linpred() call left is the covariance
  # check, and the count no longer depends on the number of
  # coefficients at all. The next test pins the same number on a model
  # with 110 random coefficients.
  o <- sp_curve_fit(k = 10)
  g <- data.frame(x = seq(0, 1, length.out = 20))
  cv <- frm_curve(o$fit, newdata = g, simultaneous = FALSE)
  expect_equal(attr(cv, "check")$n_predict, 1L)
  expect_equal(ncol(attr(cv, "Sigma")), nrow(g))
  expect_equal(nrow(attr(cv, "Sigma")), nrow(g))
})

test_that("a grouping block costs nothing, because there is no probe", {
  # This test used to pin a probe count: the design was rebuilt by unit
  # perturbation and coefficients were screened in chunks of 24, so a
  # 40-level grouping block cost one call per chunk rather than one per
  # level. frm_lp_basis() returns the design core already had, so there
  # is no perturbation, no chunk and no probe, and the only frm_linpred()
  # call left is the covariance check. What is worth pinning now is that
  # the count does not respond to the block at all.
  set.seed(3)
  n <- 400
  d <- data.frame(x = sort(stats::runif(n)),
                  g = factor(rep(1:40, length.out = n)))
  d$y <- 2 * sin(pi * d$x) + stats::rnorm(40, 0, 0.5)[d$g] +
    stats::rnorm(n, 0, 0.4)
  fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 9) + (1 | g)),
                     family = stats::gaussian(), data = d)
  g <- data.frame(x = seq(0, 1, length.out = 20),
                  g = factor(1, levels = levels(d$g)))
  cv <- frm_curve(fit, newdata = g, re_formula = NA, simultaneous = FALSE)
  expect_equal(attr(cv, "check")$n_predict, 1L)
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-10)
  # 40 levels in the fit and the same single call, which is the property
  # the old chunking existed to approximate
  expect_equal(length(fit$estimates$b), 47L)
})

test_that("the simultaneous band is wider than the pointwise one, and covers", {
  o <- sp_curve_fit()
  g <- data.frame(x = seq(0.02, 0.98, length.out = 50))
  cv <- frm_curve(o$fit, newdata = g, nsim = 20000, seed = 11)
  expect_gt(cv$.crit_sim[1L], cv$.crit[1L])
  expect_true(all(cv$.lower_sim <= cv$.lower_ci))
  expect_true(all(cv$.upper_sim >= cv$.upper_ci))
  # the critical value comes back with its own Monte Carlo standard
  # error, and at 20000 draws that error is small but not zero
  mcse <- attr(cv, "check")$crit_mcse
  expect_true(is.finite(mcse))
  expect_lt(mcse, 0.05)
  expect_gt(mcse, 0)
  # the simultaneous band covers the truth everywhere; the pointwise one
  # is not asked to
  truth <- 2 * sin(pi * g$x) + 0.6 * g$x
  expect_true(all(truth >= cv$.lower_sim & truth <= cv$.upper_sim))
})

test_that("the critical value is stable in the seed and moves with nsim", {
  o <- sp_curve_fit()
  g <- data.frame(x = seq(0, 1, length.out = 30))
  a <- frm_curve(o$fit, newdata = g, nsim = 50000, seed = 2)
  b <- frm_curve(o$fit, newdata = g, nsim = 50000, seed = 2)
  expect_identical(a$.crit_sim[1L], b$.crit_sim[1L])
  cc <- frm_curve(o$fit, newdata = g, nsim = 50000, seed = 99)
  # two seeds differ by simulation noise alone. The difference of two
  # independent quantile estimates has standard
  # deviation sqrt(2) times the reported one, so this is a four-sigma
  # bound rather than an eight-sigma one
  expect_lt(abs(a$.crit_sim[1L] - cc$.crit_sim[1L]),
            8 * attr(a, "check")$crit_mcse)
})

test_that("transform returns the band through the link inverse", {
  set.seed(5)
  d <- data.frame(x = sort(stats::runif(300)))
  d$y <- stats::rpois(300, exp(0.5 + sin(2 * pi * d$x)))
  fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 8)),
                     family = stats::poisson(), data = d)
  g <- data.frame(x = seq(0, 1, length.out = 20))
  lk <- frm_curve(fit, newdata = g, nsim = 2000, seed = 1)
  rs <- frm_curve(fit, newdata = g, nsim = 2000, seed = 1, transform = TRUE)
  expect_equal(rs$.estimate, exp(lk$.estimate), tolerance = 1e-10)
  expect_equal(rs$.lower_sim, exp(lk$.lower_sim), tolerance = 1e-10)
  expect_true(all(rs$.lower_ci > 0))
  # a transformed band is no longer symmetric, so no standard error is
  # offered on that scale
  expect_null(rs$.se)
})

test_that("the refusals name what is wrong", {
  o <- sp_curve_fit()
  g <- data.frame(x = seq(0, 1, length.out = 10))
  expect_error(frm_curve(o$d, newdata = g), "must be a frmtmb fit")
  expect_error(frm_curve(o$fit, newdata = data.frame()),
               "at least one row")
  expect_error(frm_curve(o$fit, newdata = g, level = 1.5),
               "strictly between 0 and 1")
  expect_error(frm_curve(o$fit, newdata = g, simultaneous = NA),
               "must be TRUE or FALSE")
  expect_error(frm_curve(o$fit, newdata = g[1, , drop = FALSE]),
               "at least two points")
  expect_error(frm_curve(o$fit, newdata = g, nsim = 0), "whole number")
  # a tolerance no covariance could meet refuses rather than returning
  expect_error(frm_curve(o$fit, newdata = g, tol = 0),
               "disagrees with frm_linpred")
})

test_that("a fit with no random-effect block works, through cov.fixed", {
  # This used to refuse, and the refusal was an artifact of this package
  # recomputing the joint precision itself: a model with no random
  # effects has none. Reading core's cached covariance instead picks up
  # the fallback get_joint_cov() already has (cov.fixed), so the curve
  # is simply the fixed-effect delta method and comes out exact.
  set.seed(2)
  d <- data.frame(x = stats::rnorm(80))
  d$y <- stats::rnorm(80, 1 + 2 * d$x, 0.5)
  fit <- frmtmb::frm(frmtmb::bf(y ~ x), family = stats::gaussian(), data = d)
  g <- data.frame(x = seq(-2, 2, length.out = 9))
  cv <- frm_curve(fit, newdata = g, simultaneous = FALSE)
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-10)
  p <- frm_linpred(fit, newdata = g, type = "link", se.fit = TRUE,
                      re_formula = NA)
  expect_equal(cv$.se, as.numeric(p$se.fit), tolerance = 1e-14)
  expect_true(all(cv$.se > 0))
})

test_that("the covariance is core's cached one, not a second sdreport", {
  # The cache is what makes the joint-precision solve free on the second
  # and every later call, and it is what keeps autoscaling correct. If
  # this ever stops being populated the fallback still works, so the
  # assertion is that the cache IS the source rather than that the
  # answer is right.
  o <- sp_curve_fit()
  expect_null(o$fit$cache$Vjoint)
  cv <- frm_curve(o$fit, newdata = data.frame(x = seq(0, 1, length.out = 12)),
                  simultaneous = FALSE)
  expect_false(is.null(o$fit$cache$Vjoint))
  expect_false(is.null(o$fit$cache$Vjoint$V))
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-10)
})

test_that("a t2 tensor smooth assembles like any other", {
  set.seed(4)
  n <- 400
  d <- data.frame(x = stats::runif(n), z = stats::runif(n))
  d$y <- sin(pi * d$x) * d$z + stats::rnorm(n, 0, 0.3)
  fit <- frmtmb::frm(frmtmb::bf(y ~ t2(x, z, k = 4)),
                     family = stats::gaussian(), data = d)
  g <- data.frame(x = seq(0.05, 0.95, length.out = 12), z = 0.5)
  cv <- frm_curve(fit, newdata = g, nsim = 2000, seed = 1)
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-8)
  expect_true(all(is.finite(cv$.se)))
})

test_that("a reduced-rank block is caught by the check, not by the probe", {
  # This is the row the covariance check exists for. An rr block's
  # loadings live in theta, so eta is linear in b at fixed theta and the
  # linearity probe passes; the perturbation then misses the derivative
  # with respect to the loadings, which core's delta method carries as
  # rr_jacobians(). The population prediction drops the block entirely
  # and is fine; the conditional one is refused.
  set.seed(11)
  n_site <- 30
  n_sp <- 5
  d <- expand.grid(site = factor(seq_len(n_site)), sp = factor(seq_len(n_sp)))
  d$x <- stats::runif(nrow(d))
  L <- matrix(stats::rnorm(n_sp * 2, 0, 0.7), n_sp, 2)
  fm <- matrix(stats::rnorm(n_site * 2), 2, n_site)
  d$y <- sin(pi * d$x) +
    as.vector(t(L %*% fm))[(as.integer(d$sp) - 1L) * n_site +
                             as.integer(d$site)] +
    stats::rnorm(nrow(d), 0, 0.4)
  fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 6) + rr(0 + sp | site, d = 2)),
                     family = stats::gaussian(), data = d)
  g <- data.frame(x = seq(0, 1, length.out = 15),
                  sp = factor(1, levels = levels(d$sp)),
                  site = factor(1, levels = levels(d$site)))
  cv <- frm_curve(fit, newdata = g, re_formula = NA, simultaneous = FALSE)
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-8)
  # re_formula = NULL used to be REFUSED here, and not by the linearity
  # probe: a reduced-rank block's loadings live in theta, so eta is
  # linear in b at fixed theta and the probe passed, while the design
  # the perturbation could build was missing the derivative with
  # respect to the loadings and the standard errors came out 27 percent
  # away from frm_linpred(se.fit = TRUE)'s. frm_lp_basis() carries the
  # loading columns through rr_jacobians(), so it now works.
  cvn <- frm_curve(fit, newdata = g, re_formula = NULL, simultaneous = FALSE)
  expect_lt(attr(cvn, "check")$cov_rel_error, 1e-8)
  pn <- frm_linpred(fit, newdata = g, type = "link", re_formula = NULL,
                       se.fit = TRUE)
  expect_equal(cvn$.se, as.numeric(pn$se.fit), tolerance = 1e-10)
})

test_that("a curve on a dpar other than mu finds its coefficients", {
  # Everything but the location parameter keeps its coefficients in
  # betad rather than in beta, so a design rebuilt from beta and b alone
  # would miss the fixed part of the curve entirely. It does not fail
  # quietly when that happens: the covariance check refuses. This is the
  # regression test for that.
  set.seed(8)
  d <- data.frame(x = sort(stats::runif(400)))
  d$y <- stats::rnorm(400, 1 + 0.5 * d$x, exp(-1 + sin(pi * d$x)))
  fit <- frmtmb::frm(frmtmb::bf(y ~ x, sigma ~ s(x, k = 6)),
                     family = stats::gaussian(), data = d)
  g <- data.frame(x = seq(0.05, 0.95, length.out = 20))
  cv <- frm_curve(fit, newdata = g, dpar = "sigma", nsim = 2000, seed = 1)
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-10)
  expect_true(all(cv$.se > 0))
  p <- frm_linpred(fit, newdata = g, type = "link", dpar = "sigma",
                      se.fit = TRUE, re_formula = NA)
  expect_equal(cv$.estimate, as.numeric(p$fit), tolerance = 1e-12)
  expect_equal(cv$.se, as.numeric(p$se.fit), tolerance = 1e-10)
  # and the derivative of that same dpar curve
  d1 <- frm_curve_deriv(fit, var = "x", order = 1, newdata = g,
                        dpar = "sigma", simultaneous = FALSE)
  expect_true(all(is.finite(d1$.se)))
})

test_that("a dpar held fixed does not become a design column", {
  # A fixed dpar moves the prediction when it is perturbed and has no
  # row in the joint precision, so including it would make the design
  # wider than the covariance. It is excluded before the probe runs.
  set.seed(12)
  d <- data.frame(x = sort(stats::runif(300)))
  d$y <- stats::rnorm(300, sin(pi * d$x), 0.5)
  fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 6), sigma = 0.5),
                     family = stats::gaussian(), data = d)
  cv <- frm_curve(fit, newdata = data.frame(x = seq(0, 1, length.out = 12)),
                  simultaneous = FALSE)
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-10)
})

test_that("a factor-smooth model costs the documented number of calls", {
  # The @section Cost: of frm_curve() quotes this model by name, so the
  # number in the documentation is pinned to the model it came from.
  #
  # `re_formula = NA` KEEPS every smooth from the development version
  # after 0.64.0, a factor smooth included, so the grid has to name a
  # level of `subject` and the curve is that level's. The second fit,
  # with no population smooth, used to give a CONSTANT here because the
  # fs term was dropped; it now gives that subject's own curve, and two
  # subjects give different curves. The call count, which is what this
  # block is for, is unchanged either way.
  set.seed(4)
  n_sub <- 20
  n_rep <- 12
  n_t <- 30
  peak <- function(t, h, s) h * exp(-0.5 * ((t - s) / 0.16)^2)
  sub <- rep(seq_len(n_sub), each = n_rep * n_t)
  d <- data.frame(subject = factor(sub),
                  t = rep(seq(0, 1, length.out = n_t),
                          times = n_sub * n_rep))
  h <- stats::rnorm(n_sub, 1, 0.12)
  sv <- stats::rnorm(n_sub, 0.5, 0.04)
  d$v <- peak(d$t, h[sub], sv[sub]) + stats::rnorm(nrow(d), 0, 0.06)

  # suppressWarnings: this fit reaches a maximum absolute gradient of
  # about 1.5e-3 on some platforms and not others. What this test
  # measures is the number of frm_linpred() calls and the covariance
  # identity below, neither of which that touches.
  fit <- suppressWarnings(frmtmb::frm(
    frmtmb::bf(v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)),
    family = stats::gaussian(), data = d))
  expect_equal(length(fit$estimates$b), 110L)
  g <- data.frame(t = seq(0, 1, length.out = 80),
                  subject = factor(levels(d$subject)[1],
                                   levels = levels(d$subject)))
  cv <- frm_curve(fit, newdata = g, re_formula = NA, simultaneous = FALSE)
  # 32 frm_linpred() calls before the seam; one now, on a model with 110
  # random coefficients, which is the point of the seam
  expect_equal(attr(cv, "check")$n_predict, 1L)
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-10)
  # a real bell, not a constant
  expect_gt(stats::sd(cv$.estimate), 0.2)
  # the grouping column is required, and says so by name
  expect_error(frm_curve(fit, newdata = g[, "t", drop = FALSE],
                         re_formula = NA, simultaneous = FALSE),
               "needs the grouping column `subject`")

  # the same fs term with NO population smooth: the curve is that
  # subject's own, and a different subject gives a different one
  fit2 <- frmtmb::frm(
    frmtmb::bf(v ~ s(t, subject, bs = "fs", k = 5)),
    family = stats::gaussian(), data = d)
  cv2 <- frm_curve(fit2, newdata = g, re_formula = NA, simultaneous = FALSE)
  expect_equal(attr(cv2, "check")$n_predict, attr(cv, "check")$n_predict)
  expect_lt(attr(cv2, "check")$cov_rel_error, 1e-10)
  expect_gt(stats::sd(cv2$.estimate), 0.2)
  g2 <- transform(g, subject = factor(levels(d$subject)[2],
                                      levels = levels(d$subject)))
  cv3 <- frm_curve(fit2, newdata = g2, re_formula = NA, simultaneous = FALSE)
  # relative to the size of the curves themselves, not an absolute gap
  expect_gt(max(abs(cv3$.estimate - cv2$.estimate)) /
              stats::sd(cv2$.estimate), 0.1)
})

test_that("an autoscaled fit works, because the covariance is core's", {
  # An autoscaled fit carries par_units and its covariance comes from
  # autoscale_sdreport(), which reparameterizes. A fresh
  # RTMB::sdreport() would not, and this package used to call one, so
  # this path would have refused. Reading core's cache makes it correct
  # rather than merely caught.
  set.seed(7)
  n <- 300
  d <- data.frame(x = sort(stats::runif(n)), z = stats::rnorm(n, 0, 1e6))
  d$y <- 2 * sin(pi * d$x) + 1e-6 * d$z + stats::rnorm(n, 0, 0.4)
  fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 10) + z),
                     family = stats::gaussian(), data = d,
                     control = frmtmb::frmtmb_control(autoscale = TRUE))
  expect_false(is.null(fit$par_units))
  expect_false(all(fit$par_units == 1))
  g <- data.frame(x = seq(0, 1, length.out = 20), z = 0)
  cv <- frm_curve(fit, newdata = g, simultaneous = FALSE)
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-10)
  p <- frm_linpred(fit, newdata = g, type = "link", se.fit = TRUE,
                      re_formula = NA)
  expect_equal(cv$.se, as.numeric(p$se.fit), tolerance = 1e-13)
})

test_that("the covariance is core's, at the rows core says it is at", {
  # The ordering this test used to guard is gone with the machinery it
  # guarded: there is no cache read to warm and no fallback sdreport()
  # to go round autoscale_sdreport(). What is left to check is that the
  # covariance frm_lp_basis() hands over is the one it says it is: the
  # joint covariance of everything the fit estimates, subset to exactly
  # the rows the design's columns sit at.
  o <- sp_curve_fit()
  nd <- data.frame(x = c(0.2, 0.8))
  lb <- frmtmb::frm_lp_basis(o$fit, newdata = nd, re_formula = NA)
  jc <- frmtmb::frm_joint_cov(o$fit)
  expect_equal(lb$V, jc$V[lb$coef_pos, lb$coef_pos, drop = FALSE])
  expect_identical(lb$coef_names, jc$labels[lb$coef_pos])
  # and the b rows a smooth needs are in it, which is what vcov() cannot
  # reach
  expect_true(any(jc$names[lb$coef_pos] == "b"))

  cv <- frm_curve(o$fit, newdata = nd, simultaneous = FALSE)
  expect_equal(attr(cv, "Sigma"),
               unname(lb$A %*% lb$V %*% t(lb$A)), tolerance = 1e-12)
})

test_that("a band past an exact gp()'s positions draws from its kriging", {
  skip_on_cran()
  # The construction of dev/gpby-crit.R: y ~ fac + gp(x), 60 points on
  # [0, 6], a grid on [7, 12] past them. The band used to draw from
  # A V A' and divide by a standard error that carried the kriging
  # variance, so it was too narrow exactly where it should widen: its
  # critical value was 2.02211 against 2.33046 for the same draw
  # standardized by its own scale (dev/reviews/2026-09-08-diffcurve.md,
  # R9, which put the shortfall at 17 percent on its data).
  set.seed(1)
  n <- 60
  d <- data.frame(x = stats::runif(n, 0, 6),
                  fac = factor(rep(c("A", "B"), length.out = n)))
  d$y <- sin(d$x) + ifelse(d$fac == "B", 0.5, 0) + stats::rnorm(n, 0, 0.2)
  fit <- frmtmb::frm(frmtmb::bf(y ~ fac + gp(x)),
                     family = stats::gaussian(), data = d)
  nd <- data.frame(x = seq(7, 12, length.out = 25),
                   fac = factor("A", levels = levels(d$fac)))
  cv <- frm_curve(fit, newdata = nd, nsim = 20000, seed = 1)
  lb <- frmtmb::frm_lp_basis(fit, newdata = nd, re_formula = NA)
  expect_gt(min(lb$extra_var / cv$.se^2), 0.01)
  # the grid covariance is the whole of it, so its diagonal is .se^2
  expect_lt(max(abs(sqrt(diag(attr(cv, "Sigma"))) / cv$.se - 1)), 1e-10)
  # and the critical value is at least the bound that only fixes the
  # scale, which a residual decorrelating faster than the mean raises
  A <- as.matrix(lb$A)
  S0 <- A %*% lb$V %*% t(A)
  bound <- sp_sim_crit(S0, sqrt(diag(S0)), 20000, 0.95, seed = 1)
  expect_gt(cv$.crit_sim[1], bound$crit - 4 * bound$mcse)
})

test_that("a curve along a direction the fit lost gets NaN and one warning", {
  # y ~ a + b with a and b both ~ 1 + x: only a + b is determined, so a
  # curve of `a` alone has no standard error. The covariance check used
  # to stop with "this package is reading the seam wrongly" instead.
  set.seed(955)
  dd <- data.frame(x = stats::rnorm(60))
  dd$y <- 3 + 0.5 * dd$x + stats::rnorm(60, 0, 0.4)
  fit <- suppressWarnings(frmtmb::frm(
    frmtmb::bf(y ~ a + b, a ~ 1 + x, b ~ 1 + x, nl = TRUE), data = dd))
  nd <- data.frame(x = seq(-2, 2, length.out = 9))
  lost <- function(expr) {
    n <- 0L
    v <- withCallingHandlers(expr, frmtmb_se_lost_prediction = function(w) {
      n <<- n + 1L
      invokeRestart("muffleWarning")
    })
    list(value = v, n = n)
  }
  for (sim in c(FALSE, TRUE)) {
    r <- lost(frm_curve(fit, newdata = nd, dpar = "a", simultaneous = sim,
                        nsim = 500, seed = 1))
    expect_true(all(is.nan(r$value$.se)))
    expect_identical(r$n, 1L)
    r <- lost(frm_curve_deriv(fit, var = "x", newdata = nd, dpar = "a",
                              simultaneous = sim, nsim = 500, seed = 1))
    expect_true(all(is.nan(r$value$.se)))
    expect_identical(r$n, 1L)
  }
  # the mean, which the data determine, keeps its band and says nothing
  expect_no_warning(cv <- frm_curve(fit, newdata = nd, simultaneous = TRUE,
                                    nsim = 500, seed = 1))
  expect_true(all(is.finite(cv$.se)))
})
