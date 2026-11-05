# frm_linpred(type = "link"), frm_lp_basis(), frm_curve() and its inference
# on a predictor with several layers: a cs() term's thresholds and a family
# with several location predictors. Pair p013 of the equivalence survey, seed
# 13: these read the shared eta, with brms's posterior_linpred() and frmtmb's
# own fitted(scale = "linear") at eta + cs_k. Seen to fail on 0.78.0
# (rellib-r16): frm_linpred() had length n and frm_curve() n rows
# (dev/predread-cscurve.R). The brms half is in
# test-linpred-layers-brms.R.

ll_data <- function() {
  set.seed(13)
  n <- 240
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- cut(d$x - 0.5 * d$z + stats::rlogis(n),
             c(-Inf, -0.8, 0.2, 1.2, Inf), labels = FALSE)
  d
}
ll_fit <- function() {
  frm(bf(y ~ x + cs(z)), family = sratio(), data = ll_data())
}
ll_grid <- function() data.frame(x = c(-1, 0, 1, 2), z = c(0.5, 0.5, -1, 0))

test_that("frm_linpred(link) is eta + cs_k per threshold, as a matrix", {
  fit <- ll_fit()
  nd <- ll_grid()
  m <- frm_linpred(fit, newdata = nd, type = "link")
  shared <- frm_linpred(fit, newdata = nd, type = "link", dpar = "mu")
  expect_null(dim(shared))
  bx <- fixef(fit)["x", "Estimate"]
  bcs <- fixef(fit)[c("z[1]", "z[2]", "z[3]"), "Estimate"]
  expect_identical(dim(m), c(4L, 3L))
  expect_identical(colnames(m), paste0("eta", 1:3))
  expect_equal(as.numeric(shared), bx * nd$x)
  expect_equal(unname(m), unname(outer(bx * nd$x, rep(1, 3)) +
                                   outer(nd$z, unname(bcs))))
  # in sample too, on the rows of the data
  s <- frm_linpred(fit, type = "link")
  expect_identical(dim(s), c(nrow(ll_data()), 3L))
  expect_equal(unname(s[1L, ]), unname(bx * ll_data()$x[1L] +
                                         ll_data()$z[1L] * bcs))
})

test_that("se.fit has the matrix's shape and reads the cs() coefficients", {
  fit <- ll_fit()
  nd <- ll_grid()
  p <- frm_linpred(fit, newdata = nd, type = "link", se.fit = TRUE)
  expect_identical(dim(p$se.fit), dim(p$fit))
  # by hand: var(b_x x + bcs_k z), from the joint covariance
  jc <- frm_joint_cov(fit)
  at <- which(jc$labels == "beta.x")
  ck <- which(grepl("^bcs", jc$labels))
  expect_length(ck, 3L)
  for (k in 1:3) {
    a <- cbind(nd$x, nd$z)
    i <- c(at, ck[k])
    v <- rowSums((a %*% jc$V[i, i]) * a)
    expect_equal(unname(p$se.fit[, k]), sqrt(v))
  }
  # fitted(scale = "linear") reads the same builder, and its finite-difference
  # route of 0.78.0 is one more reading of the same number
  fi <- fitted(fit, newdata = nd, scale = "linear")
  expect_equal(unname(fi[, "Estimate", ]), unname(p$fit))
  expect_equal(unname(fi[, "Est.Error", ]), unname(p$se.fit))
  fd <- fit_fd_se(fit, function(f) {
    frm_linpred(f, newdata = nd, type = "link")
  })
  expect_lt(max(abs(fd / p$se.fit - 1)), 1e-4)
})

test_that("a model with no second predictor keeps its vector and its frame", {
  set.seed(2)
  d <- ll_data()
  fit <- frm(bf(y ~ x + z), family = sratio(), data = d)
  nd <- ll_grid()
  v <- frm_linpred(fit, newdata = nd, type = "link")
  expect_null(dim(v))
  expect_false(".eta" %in% names(frm_curve(fit, newdata = nd,
                                           simultaneous = FALSE)))
  expect_null(frm_lp_basis(fit, newdata = nd)$layers)
})

test_that("dpar = 'mu' is the shared part on every reader, as in brms", {
  fit <- ll_fit()
  nd <- ll_grid()
  shared <- as.numeric(frm_linpred(fit, newdata = nd, dpar = "mu"))
  expect_equal(frm_curve(fit, newdata = nd, dpar = "mu",
                         simultaneous = FALSE)$.estimate, shared)
  expect_equal(as.numeric(frm_lp_basis(fit, newdata = nd, dpar = "mu")$eta),
               shared)
  # brms 2.23.0, fitted(dpar = "mu", scale = "linear"): 3 x 4, no threshold
  # layers (dev/predread-cs-layout-brms.R); on 0.78.0 frmtmb gave them
  fi <- fitted(fit, newdata = nd, dpar = "mu", scale = "linear")
  expect_length(dim(fi), 2L)
  expect_equal(unname(fi[, "Estimate"]), shared)
})

test_that("dpar = 'mu' on the response scale is the family's link inverse", {
  # brms 2.23.0's posterior_epred(dpar = "mu") and
  # posterior_linpred(transform = TRUE); on 0.78.0 these read the shared part
  # unchanged, because the location's own link slot is the identity
  fit <- ll_fit()
  nd <- ll_grid()
  sh <- as.numeric(frm_linpred(fit, newdata = nd, dpar = "mu"))
  r <- frm_linpred(fit, newdata = nd, dpar = "mu", type = "response")
  expect_equal(as.numeric(r), stats::plogis(sh))
  expect_equal(unname(fitted(fit, newdata = nd, dpar = "mu")[, "Estimate"]),
               stats::plogis(sh))
  p <- frm_linpred(fit, newdata = nd, dpar = "mu", type = "response",
                   se.fit = TRUE)
  b <- frm_linpred(fit, newdata = nd, dpar = "mu", se.fit = TRUE)
  expect_equal(as.numeric(p$se.fit),
               as.numeric(b$se.fit) * stats::dlogis(sh))
  cv <- frm_curve(fit, newdata = nd, dpar = "mu", transform = TRUE,
                  simultaneous = FALSE)
  expect_equal(cv$.estimate, stats::plogis(sh))
  ce <- conditional_effects(fit, "x", dpar = "mu", resolution = 3L)[[1L]]
  expect_equal(ce$estimate__,
               stats::plogis(fixef(fit)["x", "Estimate"] * ce$x))
  # a probit family reads its own link
  fp <- frm(bf(y ~ x), family = cumulative("probit"), data = ll_data())
  sp <- as.numeric(frm_linpred(fp, newdata = nd, dpar = "mu"))
  expect_equal(as.numeric(frm_linpred(fp, newdata = nd, dpar = "mu",
                                      type = "response")), stats::pnorm(sp))
  # the density still reads the predictor: the probabilities are unchanged
  P <- frm_linpred(fit, newdata = nd, type = "response")
  expect_equal(unname(rowSums(P)), rep(1, nrow(nd)))
})

test_that("frm_curve() returns one curve per threshold, jointly covered", {
  fit <- ll_fit()
  nd <- ll_grid()
  cv <- frm_curve(fit, newdata = nd, nsim = 2000, seed = 3)
  m <- frm_linpred(fit, newdata = nd, type = "link")
  expect_identical(nrow(cv), 12L)
  expect_identical(levels(cv$.eta), paste0("eta", 1:3))
  expect_identical(as.character(cv$.eta), rep(paste0("eta", 1:3), each = 4L))
  expect_equal(cv$x, rep(nd$x, 3L))
  expect_equal(cv$.estimate, as.numeric(m))
  expect_equal(cv$.se, as.numeric(frm_linpred(fit, newdata = nd,
                                              se.fit = TRUE)$se.fit))
  S <- attr(cv, "Sigma")
  expect_identical(dim(S), c(12L, 12L))
  # a row's two thresholds share b_x, so their covariance is its variance
  # plus the z part they do not share; at z = 0 (row 4) it is that variance
  jc <- frm_joint_cov(fit)
  vx <- jc$V[jc$labels == "beta.x", jc$labels == "beta.x"]
  expect_equal(S[4L, 8L], nd$x[4L]^2 * vx)
  # the check ran over the layers
  expect_true(is.finite(attr(cv, "check")$cov_rel_error))
  # a band over twelve correlated points is wider than one point's
  expect_gt(cv$.crit_sim[1L], cv$.crit[1L])
})

test_that("a difference of two grids is taken layer by layer", {
  fit <- ll_fit()
  nd <- ll_grid()
  ct <- transform(nd, x = x + 1)
  dc <- frm_curve(fit, newdata = nd, contrast = ct, simultaneous = FALSE)
  m <- frm_linpred(fit, newdata = nd) - frm_linpred(fit, newdata = ct)
  expect_equal(dc$.estimate, as.numeric(m))
  jc <- frm_joint_cov(fit)
  vx <- jc$V[jc$labels == "beta.x", jc$labels == "beta.x"]
  # only b_x separates the grids, on every layer
  expect_equal(dc$.se, rep(sqrt(vx), 12L))
})

test_that("the derivative of a cs() variable differs by threshold", {
  fit <- ll_fit()
  nd <- ll_grid()
  dz <- frm_curve_deriv(fit, var = "z", newdata = nd, simultaneous = FALSE)
  fz <- fixef(fit)
  expect_identical(nrow(dz), 12L)
  expect_equal(dz$.estimate, rep(unname(fz[paste0("z[", 1:3, "]"),
                                           "Estimate"]), each = 4L),
               tolerance = 1e-6)
  expect_equal(dz$.se, rep(unname(fz[paste0("z[", 1:3, "]"), "Est.Error"]),
                           each = 4L), tolerance = 1e-6)
  dx <- frm_curve_deriv(fit, var = "x", newdata = nd, simultaneous = FALSE)
  expect_equal(dx$.estimate, rep(unname(fz["x", "Estimate"]), 12L),
               tolerance = 1e-6)
  # the second order of a linear predictor is zero on every layer
  d2 <- frm_curve_deriv(fit, var = "x", order = 2, newdata = nd,
                        simultaneous = FALSE)
  expect_identical(nrow(d2), 12L)
  expect_lt(max(abs(d2$.estimate)), 1e-4 * abs(fz["x", "Estimate"]))
})

test_that("a feature is searched on each curve", {
  fit <- ll_fit()
  # no grid point at the root: the scan needs a change of sign
  grid <- data.frame(x = seq(-3, 3, length.out = 14), z = 1)
  ft <- frm_curve_feature(fit, var = "x", type = "crossing", at = 0,
                          newdata = grid)
  fz <- fixef(fit)
  # eta + cs_k = b_x x + bcs_k * 1 crosses zero at -bcs_k / b_x
  root <- -fz[paste0("z[", 1:3, "]"), "Estimate"] / fz["x", "Estimate"]
  expect_equal(ft$.estimate, unname(root), tolerance = 1e-6)
  expect_identical(as.character(ft$.eta), paste0("eta", 1:3))
  expect_true(all(is.finite(ft$.se)))
  # dpar = "mu" is one curve, crossing at x = 0
  ft0 <- frm_curve_feature(fit, var = "x", type = "crossing", at = 0,
                           newdata = grid, dpar = "mu")
  expect_false(".eta" %in% names(ft0))
  expect_equal(ft0$.estimate, 0, tolerance = 1e-6)
})

test_that("a curve refuses a cov = FALSE term with its reason", {
  set.seed(6)
  d <- data.frame(t = 1:60, x = stats::rnorm(60))
  d$y <- 1 + 0.5 * d$x + stats::arima.sim(list(ar = 0.5), 60)
  d$y <- as.numeric(d$y)
  fit <- frm(bf(y ~ x + ar(time = t, cov = FALSE)), data = d,
             family = gaussian())
  nd <- data.frame(x = c(-1, 1), t = 61:62, y = c(0, 1))
  expect_error(frm_curve(fit, newdata = nd, simultaneous = FALSE),
               "reads the residual of the row before")
  # the one-step mean is frm_linpred()'s; sigma has a curve of its own
  expect_length(frm_linpred(fit, newdata = nd), 2L)
})

test_that("a new level's draw is in every layer, and covaries across them", {
  d <- ll_data()
  d$g <- factor(rep(1:12, length.out = nrow(d)))
  fit <- allow_warnings(
    frm(bf(y ~ x + cs(z) + (1 | g)), family = sratio(), data = d),
    c("variance", "boundary", "gradient"))
  nd <- data.frame(x = c(0.1, 1.1), z = 0.5, g = factor("new"))
  cv <- frm_curve(fit, newdata = nd, re_formula = NULL,
                  allow_new_levels = TRUE, simultaneous = FALSE)
  expect_identical(nrow(cv), 6L)
  S <- attr(cv, "Sigma")
  lb <- frm_lp_basis(fit, newdata = nd, re_formula = NULL,
                     allow_new_levels = TRUE, extra_cov = TRUE)
  E <- as.matrix(lb$extra_cov)
  expect_identical(dim(E), c(6L, 6L))
  # one draw per unseen level, whichever threshold reads it
  expect_equal(E[1L, 4L], E[1L, 1L])
  expect_gt(E[1L, 1L], 0)
  # the slope of a random intercept's draw is 0, so the derivative's standard
  # error is the coefficient's alone
  dx <- frm_curve_deriv(fit, var = "x", newdata = nd, re_formula = NULL,
                        allow_new_levels = TRUE, simultaneous = FALSE)
  expect_equal(dx$.se, rep(unname(fixef(fit)["x", "Est.Error"]), 6L),
               tolerance = 1e-5)
})

test_that("emmeans averages the layers, and rep.meas names them", {
  skip_if_not_installed("emmeans")
  fit <- ll_fit()
  at <- list(x = c(-1, 1), z = 0.5)
  by_k <- as.data.frame(emmeans::emmeans(fit, ~ x | rep.meas, at = at))
  expect_identical(nrow(by_k), 6L)
  m <- frm_linpred(fit, newdata = data.frame(x = c(-1, 1), z = 0.5))
  expect_equal(by_k$emmean, as.numeric(m))
  p <- frm_linpred(fit, newdata = data.frame(x = c(-1, 1), z = 0.5),
                   se.fit = TRUE)
  expect_equal(by_k$SE, as.numeric(p$se.fit))
  av <- as.data.frame(emmeans::emmeans(fit, ~ x, at = at))
  expect_equal(av$emmean, unname(rowMeans(m)))
  # dpar = "mu": the shared part, no rep.meas, on the design route
  sh <- as.data.frame(emmeans::emmeans(fit, ~ x, at = at, dpar = "mu"))
  expect_equal(sh$emmean, c(-1, 1) * fixef(fit)["x", "Estimate"])
})

test_that("marginaleffects' get_predict keys the layers in `group`", {
  skip_if_not_installed("marginaleffects")
  fit <- ll_fit()
  nd <- ll_grid()
  g <- marginaleffects::get_predict(fit, newdata = nd, type = "link")
  expect_identical(g$group, rep(paste0("eta", 1:3), each = 4L))
  expect_equal(g$estimate, as.numeric(frm_linpred(fit, newdata = nd)))
})

test_that("insight still reads one number per row", {
  skip_if_not_installed("insight")
  fit <- ll_fit()
  nd <- ll_grid()
  gp <- insight::get_predicted(fit, data = nd, predict = "link", ci = NULL)
  expect_length(gp, nrow(nd))
  expect_equal(as.numeric(gp), as.numeric(frm_linpred(fit, newdata = nd,
                                                      dpar = "mu")))
})

test_that("the predictive readers keep the shared part beside .cs", {
  # the densities read tau_k - (eta + cs_k): a stacked eta here would count
  # the cs() part twice, or fail on a matrix mu. Measured against the
  # category probabilities, which never read frm_linpred(link).
  fit <- ll_fit()
  nd <- ll_grid()
  P <- frm_linpred(fit, newdata = nd, type = "response")
  set.seed(4)
  dr <- predict(fit, newdata = nd, summary = FALSE, ndraws = 4000L,
                propagate_error = FALSE)
  for (i in seq_len(nrow(nd))) {
    obs <- tabulate(dr[, i], 4L) / 4000
    se <- sqrt(P[i, ] * (1 - P[i, ]) / 4000)
    expect_lt(max(abs(obs - P[i, ]) / pmax(se, 1e-3)), 4.5)
  }
  thr <- ord_thres_linpred(fit, newdata = nd)
  tau <- fixef(fit)[paste0("Intercept[", 1:3, "]"), "Estimate"]
  m <- frm_linpred(fit, newdata = nd)
  expect_equal(unname(thr), unname(outer(rep(1, 4L), unname(tau)) - m))
})

test_that("several location predictors stack, and cs() adds its thresholds", {
  set.seed(14)
  d <- ll_data()
  fit <- allow_warnings(
    frm(bf(y ~ 1, mu1 ~ x + cs(z), mu2 ~ x),
        family = mixture(cumulative(), cumulative()), data = d,
        control = frmtmb_control(check_se = "ignore")),
    c("Category specific effects", "did not report convergence",
      "Large maximum absolute gradient"))
  nd <- ll_grid()
  m <- frm_linpred(fit, newdata = nd, type = "link")
  expect_identical(colnames(m), paste0("eta", 1:4))
  mu1 <- as.numeric(frm_linpred(fit, newdata = nd, dpar = "mu1"))
  mu2 <- as.numeric(frm_linpred(fit, newdata = nd, dpar = "mu2"))
  expect_equal(unname(m[, 4L]), mu2)
  ct <- fit$frame[["linpreds"]][["y.mu1"]][["cs"]][[1L]]
  bcs <- unname(fit$estimates[[ct$par]])
  expect_equal(unname(m[, 1:3]), unname(mu1 + outer(nd$z, bcs)))
  # brms 2.23.0 gives this fit 3 + 1 layers named eta1 to eta4
  # (dev/predread-cs-layout-brms.R)
  fi <- fitted(fit, newdata = nd, scale = "linear")
  expect_identical(dimnames(fi)[[3L]], paste0("eta", 1:4))
  expect_equal(unname(fi[, "Estimate", ]), unname(m))
  cv <- frm_curve(fit, newdata = nd, simultaneous = FALSE)
  expect_equal(cv$.estimate, as.numeric(m))

  dc <- ll_data()
  dc$y <- factor(dc$y)
  fc <- frm(bf(y ~ x), family = categorical(), data = dc)
  mc <- frm_linpred(fc, newdata = nd, type = "link")
  expect_identical(ncol(mc), 3L)
  dp1 <- fc$spec$responses[[1L]]$primary_dpars[1L]
  fc1 <- as.numeric(frm_linpred(fc, newdata = nd, dpar = dp1))
  expect_equal(unname(mc[, 1L]), fc1)
  se <- frm_linpred(fc, newdata = nd, se.fit = TRUE)$se.fit
  expect_identical(dim(se), dim(mc))
  expect_equal(frm_curve(fc, newdata = nd, simultaneous = FALSE)$.se,
               as.numeric(se))
})
