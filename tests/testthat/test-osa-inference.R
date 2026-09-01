# One-step-ahead residuals on censored and ordinal fits, and the
# joint-precision guard shared by vcov()/confint() under REML/profile.
#
# The calibration checks follow test-trunc-postfit.R: fit a correctly
# specified model, then require the residuals to be standard normal
# (KS) and to match the analytic transform they claim to be.

fit_cens_right <- function() {
  set.seed(101)
  n <- 300
  x <- stats::rnorm(n)
  ystar <- 1 + 0.8 * x + stats::rnorm(n, 0, 1.2)
  dd <- data.frame(y = pmin(ystar, 2), x = x,
                   cen = as.numeric(ystar > 2))
  list(fit = frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd),
       data = dd, cpoint = 2)
}

test_that("OSA residuals on a right-censored fit are the conditional PIT", {
  z <- fit_cens_right()
  r <- residuals(z$fit, type = "osa")
  cen <- z$data$cen != 0
  # a censored row observes an event, not a value: no residual
  expect_equal(length(r), nrow(z$data))
  expect_true(all(is.na(r[cen])))
  expect_true(all(is.finite(r[!cen])))
  expect_gt(sum(cen), 40)

  # uncensored rows are the draws that fell inside the censoring window,
  # so their PIT renormalizes on it: F(y) / F(c)
  dp <- frmtmb:::eval_dpars(z$fit)[["y"]]
  Fc <- stats::pnorm((z$cpoint - dp$mu) / dp$sigma)
  ref <- stats::qnorm(stats::pnorm((z$data$y - dp$mu) / dp$sigma) / Fc)
  expect_vector_equal(r[!cen], ref[!cen], tol = 1e-6)
  expect_gt(stats::ks.test(r[!cen], "pnorm")$p.value, 0.01)
  # ignoring the window shrinks the residuals; that was the old defect
  expect_lt(mean(stats::qnorm(
    stats::pnorm((z$data$y[!cen] - dp$mu[!cen]) / dp$sigma[!cen]))), -0.3)
})

test_that("OSA residuals survive two-sided censoring", {
  set.seed(102)
  n <- 300
  x <- stats::rnorm(n)
  ys <- 0.5 * x + stats::rnorm(n)
  lo <- -1
  hi <- 1.5
  dd <- data.frame(y = pmin(pmax(ys, lo), hi), x = x,
                   cen = ifelse(ys < lo, -1, ifelse(ys > hi, 1, 0)))
  fit <- frm(bf(y | cens(cen) ~ x) + gaussian(), data = dd)
  r <- residuals(fit, type = "osa")
  expect_true(all(is.na(r[dd$cen != 0])))
  expect_gt(stats::ks.test(r[dd$cen == 0], "pnorm")$p.value, 0.01)
})

test_that("OSA refuses the censoring layouts it cannot transform", {
  z <- fit_cens_right()
  dd <- z$data
  # row-varying censoring points: no single observation window
  dv <- dd
  i <- which(dv$cen == 1)
  dv$y[i] <- dv$y[i] + seq_along(i) * 0.01
  fv <- frm(bf(y | cens(cen) ~ x) + gaussian(), data = dv)
  expect_error(residuals(fv, type = "osa"), "type-I censoring")
  # and it names why dharma_residuals() is not the fallback
  expect_error(residuals(fv, type = "osa"), "LATENT uncensored")

  # interval censoring
  di <- dd
  di$cen[1:20] <- 2
  di$y2 <- di$y + 1
  fi <- frm(bf(y | cens(cen, y2) ~ x) + gaussian(), data = di)
  expect_error(residuals(fi, type = "osa"), "interval censoring")

  # every method that differentiates the observation hits the point mass
  expect_error(residuals(z$fit, type = "osa", osa_method = "fullGaussian"),
               "oneStepGeneric")
})

test_that("OSA residuals on ordinal fits are randomized quantile residuals", {
  set.seed(7)
  n <- 300
  x <- stats::rnorm(n)
  eta <- 0.9 * x
  thr <- c(-1, 0.3, 1.4)
  P <- rbind(stats::plogis(thr[1] - eta),
             stats::plogis(thr[2] - eta) - stats::plogis(thr[1] - eta),
             stats::plogis(thr[3] - eta) - stats::plogis(thr[2] - eta),
             1 - stats::plogis(thr[3] - eta))
  yy <- apply(P, 2, function(p) sample.int(4, 1, prob = p))
  dd <- data.frame(y = ordered(yy), x = x)

  for (fam in list(cumulative(), sratio(), cratio(), acat())) {
    fit <- frm(bf(y ~ x) + fam, data = dd)
    r <- residuals(fit, type = "osa")
    expect_true(all(is.finite(r)), info = fam$family)
    expect_gt(stats::ks.test(r, "pnorm")$p.value, 0.01)
  }

  # the taped pmf is exact: reproduce the residual analytically, using
  # oneStepPredict's own default randomization seed
  fit <- frm(bf(y ~ x) + cumulative(), data = dd)
  r <- residuals(fit, type = "osa")
  dp <- frmtmb:::eval_dpars(fit)[["y"]]
  raw <- fit$estimates[["tau_raw"]]
  tau <- cumsum(c(raw[1], exp(raw[-1])))
  cp <- rbind(stats::plogis(tau[1] - dp$mu), stats::plogis(tau[2] - dp$mu),
              stats::plogis(tau[3] - dp$mu), rep(1, n))
  Fy <- cp[cbind(yy, seq_len(n))]
  Fprev <- ifelse(yy == 1, 0, cp[cbind(pmax(yy - 1, 1), seq_len(n))])
  set.seed(123)
  ref <- stats::qnorm(Fy - stats::runif(n) * (Fy - Fprev))
  expect_vector_equal(r, ref, tol = 1e-8)
})

test_that("ordinal OSA residuals work with random effects", {
  set.seed(9)
  ng <- 25
  nt <- 12
  n <- ng * nt
  g <- factor(rep(seq_len(ng), each = nt))
  x <- stats::rnorm(n)
  eta <- 0.8 * x + stats::rnorm(ng, 0, 0.7)[g]
  P <- rbind(stats::plogis(-0.8 - eta),
             stats::plogis(0.5 - eta) - stats::plogis(-0.8 - eta),
             1 - stats::plogis(0.5 - eta))
  y <- apply(P, 2, function(p) sample.int(3, 1, prob = p))
  fit <- frm(bf(y ~ x + (1 | g)) + cumulative(),
             data = data.frame(y = ordered(y), x = x, g = g))
  r <- residuals(fit, type = "osa")
  expect_true(all(is.finite(r)))
  expect_gt(stats::ks.test(r, "pnorm")$p.value, 0.01)
})

test_that("a singular joint precision degrades like the ML branch", {
  # ML reads an already-inverted cov.fixed and reports NaN; the
  # REML/profile branch used to let solve() throw a raw LAPACK message
  Q <- matrix(c(1, 1, 1, 1), 2, 2,
              dimnames = list(c("beta", "beta"), c("beta", "beta")))
  expect_warning(V <- frmtmb:::solve_joint_precision(Q), "diagnose\\(\\)")
  expect_true(all(is.nan(V)))
  expect_equal(dim(V), c(2L, 2L))
  expect_equal(rownames(V), c("beta", "beta"))

  Qs <- Matrix::Matrix(matrix(c(1, 1, 1, 1), 2, 2), sparse = TRUE)
  expect_warning(Vs <- frmtmb:::solve_joint_precision(Qs), "diagnose\\(\\)")
  expect_true(all(is.nan(Vs)))

  # a well-conditioned precision is untouched
  Qg <- diag(c(2, 4))
  expect_silent(Vg <- frmtmb:::solve_joint_precision(Qg))
  expect_equal(unname(diag(Vg)), c(0.5, 0.25))
})
