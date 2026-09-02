# Within-group residual correlation (R-side): ar(), ma(), arma(),
# cosy(), unstr().
#
# The gate for this feature is nlme: gls() is the exact reference for a
# correlated residual with no random effect, and lme() for one with. A
# structure that cannot reproduce those log-likelihoods is not shipped,
# so the agreement tests come first and the plumbing tests after.

skip_no_nlme <- function() testthat::skip_if_not_installed("nlme")

# nlme's own optimizer defaults leave five or six digits; tightening
# them is what makes a 1e-9 comparison meaningful rather than a
# comparison of two convergence tolerances.
gls_ctl <- function() nlme::glsControl(tolerance = 1e-9, msTol = 1e-10)

ac_sim <- function(seed = 7, G = 40, K = 6, ar = 0.6, ma = NULL,
                   sd = 1.5) {
  set.seed(seed)
  d <- expand.grid(week = 1:K, subj = factor(1:G))
  e <- as.vector(vapply(seq_len(G), function(i) {
    as.vector(stats::arima.sim(list(ar = ar, ma = ma), K,
                               n.start = 200, sd = 1)) * sd
  }, numeric(K)))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + e
  d
}

# natural-scale correlation parameters of the fit's residual block
ac_pars <- function(fit) {
  ac <- fit$frame$autocor[[1L]]
  unname(autocor_natural(fit$estimates$thetaac[ac$theta_idx], ac))
}

# ------------------------------------------------------------ gls agreement

test_that("ar(p = 1) reproduces gls(corAR1) under ML and REML", {
  skip_no_nlme()
  d <- ac_sim()
  for (m in c("ML", "REML")) {
    f <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(),
             data = d, REML = (m == "REML"))
    g <- nlme::gls(y ~ x, data = d,
                   correlation = nlme::corAR1(form = ~ week | subj),
                   method = m, control = gls_ctl())
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)),
                 tolerance = 1e-9)
    expect_equal(unname(f$estimates$beta), unname(coef(g)),
                 tolerance = 1e-6)
    expect_equal(sigma(f), g$sigma, tolerance = 1e-6)
    expect_equal(ac_pars(f),
                 unname(coef(g$modelStruct$corStruct,
                             unconstrained = FALSE)),
                 tolerance = 1e-6)
  }
})

test_that("ar(p = 2) reproduces gls(corARMA(p = 2))", {
  skip_no_nlme()
  d <- ac_sim(seed = 12, ar = c(0.6, -0.3))
  for (m in c("ML", "REML")) {
    f <- frm(bf(y ~ x + ar(week, subj, p = 2, cov = TRUE)) + gaussian(),
             data = d, REML = (m == "REML"))
    g <- nlme::gls(y ~ x, data = d,
                   correlation = nlme::corARMA(form = ~ week | subj,
                                               p = 2, q = 0),
                   method = m, control = gls_ctl())
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)),
                 tolerance = 1e-8)
    expect_equal(ac_pars(f),
                 unname(coef(g$modelStruct$corStruct,
                             unconstrained = FALSE)),
                 tolerance = 1e-5)
  }
})

test_that("ma(q = 1) and arma(1, 1) reproduce gls(corARMA)", {
  skip_no_nlme()
  d3 <- ac_sim(seed = 21, ar = numeric(0), ma = 0.7)
  for (m in c("ML", "REML")) {
    f <- frm(bf(y ~ x + ma(week, subj, cov = TRUE)) + gaussian(),
             data = d3, REML = (m == "REML"))
    g <- nlme::gls(y ~ x, data = d3,
                   correlation = nlme::corARMA(0.5, form = ~ week | subj,
                                               p = 0, q = 1),
                   method = m, control = gls_ctl())
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)),
                 tolerance = 1e-9)
    expect_equal(ac_pars(f),
                 unname(coef(g$modelStruct$corStruct,
                             unconstrained = FALSE)),
                 tolerance = 1e-6)
  }
  d4 <- ac_sim(seed = 33, ar = 0.6, ma = 0.4)
  for (m in c("ML", "REML")) {
    f <- frm(bf(y ~ x + arma(week, subj, cov = TRUE)) + gaussian(),
             data = d4, REML = (m == "REML"))
    g <- nlme::gls(y ~ x, data = d4,
                   correlation = nlme::corARMA(c(0.5, 0.3),
                                               form = ~ week | subj,
                                               p = 1, q = 1),
                   method = m, control = gls_ctl())
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)),
                 tolerance = 1e-8)
    expect_equal(ac_pars(f),
                 unname(coef(g$modelStruct$corStruct,
                             unconstrained = FALSE)),
                 tolerance = 1e-5)
  }
})

test_that("cosy() and unstr() reproduce gls(corCompSymm / corSymm)", {
  skip_no_nlme()
  set.seed(5)
  G <- 40
  d <- expand.grid(week = 1:5, subj = factor(1:G))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(G, 0, 1.2)[d$subj] +
    stats::rnorm(nrow(d))
  for (m in c("ML", "REML")) {
    f <- frm(bf(y ~ x + cosy(week, subj)) + gaussian(), data = d,
             REML = (m == "REML"))
    g <- nlme::gls(y ~ x, data = d,
                   correlation = nlme::corCompSymm(form = ~ week | subj),
                   method = m, control = gls_ctl())
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)),
                 tolerance = 1e-9)
    expect_equal(ac_pars(f),
                 unname(coef(g$modelStruct$corStruct,
                             unconstrained = FALSE)),
                 tolerance = 1e-6)
  }
  for (m in c("ML", "REML")) {
    f <- suppressWarnings(
      frm(bf(y ~ x + unstr(week, subj)) + gaussian(), data = d,
          REML = (m == "REML")))
    g <- nlme::gls(y ~ x, data = d,
                   correlation = nlme::corSymm(form = ~ week | subj),
                   method = m, control = gls_ctl())
    # ten free correlations on a flat surface: the log-likelihoods pin
    # the optimum down to 1e-8, the individual correlations only to
    # about 1e-5, which is the width of the flat region both optimizers
    # stop inside
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)),
                 tolerance = 1e-7)
    expect_equal(sort(ac_pars(f)),
                 sort(unname(coef(g$modelStruct$corStruct,
                                  unconstrained = FALSE))),
                 tolerance = 1e-3)
  }
})

test_that("ragged groups keep the gap-aware lag (gls corAR1)", {
  skip_no_nlme()
  d <- ac_sim()
  set.seed(99)
  dr <- d[-sample(nrow(d), 40), ]
  for (m in c("ML", "REML")) {
    f <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(),
             data = dr, REML = (m == "REML"))
    g <- nlme::gls(y ~ x, data = dr,
                   correlation = nlme::corAR1(form = ~ week | subj),
                   method = m, control = gls_ctl())
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)),
                 tolerance = 1e-8)
    expect_equal(ac_pars(f),
                 unname(coef(g$modelStruct$corStruct,
                             unconstrained = FALSE)),
                 tolerance = 1e-5)
  }
  # more than one pattern really is exercised
  expect_gt(length(frm(bf(y ~ x + ar(week, subj, cov = TRUE)) +
                         gaussian(), data = dr,
                       dry_run = "frame")$autocor[[1L]]$patterns), 1L)
})

test_that("heterogeneous sigma matches gls(weights = varIdent)", {
  skip_no_nlme()
  d <- ac_sim(seed = 4)
  d$grp <- factor(rep(c("a", "b"), each = nrow(d) / 2))
  d$y <- d$y * ifelse(d$grp == "b", 2, 1)
  for (m in c("ML", "REML")) {
    f <- frm(bf(y ~ x + ar(week, subj, cov = TRUE), sigma ~ 0 + grp) +
               gaussian(), data = d, REML = (m == "REML"))
    g <- nlme::gls(y ~ x, data = d,
                   correlation = nlme::corAR1(form = ~ week | subj),
                   weights = nlme::varIdent(form = ~ 1 | grp),
                   method = m, control = gls_ctl())
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)),
                 tolerance = 1e-8)
    # the sigma dpar is log-linked, so exp() of its coefficients are the
    # two residual SDs varIdent reports as sigma times its ratios
    expect_equal(sort(exp(unname(f$estimates$betad))),
                 sort(unname(g$sigma *
                               coef(g$modelStruct$varStruct,
                                    unconstrained = FALSE,
                                    allCoef = TRUE))),
                 tolerance = 1e-5)
  }
})

test_that("a random intercept alongside the residual matches lme", {
  skip_no_nlme()
  d <- ac_sim(seed = 8, ar = 0.5, sd = 1)
  set.seed(8)
  d$y <- d$y + stats::rnorm(40, 0, 1)[d$subj]
  for (m in c("ML", "REML")) {
    f <- frm(bf(y ~ x + (1 | subj) + ar(week, subj, cov = TRUE)) +
               gaussian(), data = d, REML = (m == "REML"))
    g <- nlme::lme(y ~ x, random = ~ 1 | subj, data = d,
                   correlation = nlme::corAR1(form = ~ week | subj),
                   method = m,
                   control = nlme::lmeControl(msTol = 1e-10,
                                              tolerance = 1e-9))
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)),
                 tolerance = 1e-8)
    expect_equal(exp(unname(f$estimates$theta[1])),
                 as.numeric(nlme::VarCorr(g)[1, 2]), tolerance = 1e-5)
    expect_equal(sigma(f), g$sigma, tolerance = 1e-5)
    expect_equal(ac_pars(f),
                 unname(coef(g$modelStruct$corStruct,
                             unconstrained = FALSE)),
                 tolerance = 1e-5)
  }
})

test_that("the R-side ar1 is the limit of the G-side ar1 emulation", {
  # An ar1() random effect over the time factor plus a vanishing
  # independent residual approaches the R-side residual: the same
  # covariance, once the nugget is gone. A demonstration, so the
  # tolerance is loose and only the trend has to hold.
  d <- ac_sim(seed = 3, G = 25, K = 5, sd = 1.2)
  d$wk <- factor(d$week)
  ref <- as.numeric(logLik(
    frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(), data = d)))
  gaps <- vapply(c(0.2, 0.05, 0.01), function(s) {
    fit <- suppressWarnings(
      frm(bf(y ~ x + ar1(wk + 0 | subj), sigma = s) + gaussian(),
          data = d))
    abs(as.numeric(logLik(fit)) - ref)
  }, numeric(1))
  expect_true(gaps[3] < gaps[1])
  expect_lt(gaps[3], 0.5)
})

# --------------------------------------------------------------- densities

test_that("the ARMA autocorrelation function is exact", {
  # against stats::ARMAacf, which solves the same system numerically
  pac2ar <- function(p) unname(autocor_levinson(p))
  cases <- list(list(ar = 0.6, ma = numeric(0)),
                list(ar = pac2ar(c(0.5, -0.3)), ma = numeric(0)),
                list(ar = pac2ar(c(0.5, -0.3, 0.2)), ma = numeric(0)),
                list(ar = numeric(0), ma = 0.7),
                list(ar = numeric(0), ma = c(0.7, -0.4)),
                list(ar = 0.6, ma = 0.4),
                list(ar = pac2ar(c(0.5, -0.3)), ma = c(0.4, 0.2)))
  for (cs in cases) {
    ours <- autocor_arma_acf(cs$ar, cs$ma, 9)
    ref <- as.numeric(stats::ARMAacf(
      ar = if (length(cs$ar)) cs$ar else numeric(0),
      ma = if (length(cs$ma)) cs$ma else numeric(0), lag.max = 9))
    expect_equal(ours, ref, tolerance = 1e-12)
  }
})

test_that("the Levinson recursion always gives a stationary process", {
  set.seed(2)
  for (i in 1:50) {
    th <- stats::rnorm(4, 0, 3)
    phi <- autocor_levinson(autocor_pacf(th))
    roots <- polyroot(c(1, -phi))
    expect_true(all(Mod(roots) > 1 + 1e-8))
  }
})

test_that("the group densities match dmvnorm and dmvt", {
  testthat::skip_if_not_installed("mvtnorm")
  set.seed(3)
  K <- 4
  G <- 6
  ac <- list(d = K, struct = "ar", p = 1L, q = 0L,
             patterns = list(list(
               k = K, G = G, rows = seq_len(K * G),
               gather = as.vector(outer(1:K, 1:K, function(a, b) {
                 (b - 1L) * K + a
               })))))
  R <- 0.6 ^ abs(outer(1:K, 1:K, "-"))
  sg <- stats::runif(K * G, 0.5, 2)
  z <- stats::rnorm(K * G)
  refN <- 0
  refT <- 0
  for (g in seq_len(G)) {
    i <- (g - 1L) * K + seq_len(K)
    S <- diag(sg[i]) %*% R %*% diag(sg[i])
    refN <- refN + mvtnorm::dmvnorm(z[i] * sg[i], sigma = S, log = TRUE)
    refT <- refT + mvtnorm::dmvt(z[i] * sg[i], sigma = S, df = 5.5,
                                 log = TRUE)
  }
  expect_equal(autocor_loglik(z, R, ac, sum(log(sg))), refN,
               tolerance = 1e-12)
  expect_equal(autocor_loglik(z, R, ac, sum(log(sg)), nu = 5.5), refT,
               tolerance = 1e-12)
})

test_that("student() reaches the multivariate-t path", {
  d <- ac_sim(seed = 14, G = 25, K = 4, sd = 1.2)
  f <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + student(), data = d)
  expect_true(isTRUE(f$frame$autocor[[1L]]$student))
  expect_length(f$estimates$thetaac, 1L)
  expect_true(is.finite(as.numeric(logLik(f))))
  # a student fit whose nu is large sits on top of the gaussian one
  fg <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(), data = d)
  expect_gte(as.numeric(logLik(f)), as.numeric(logLik(fg)) - 1e-6)
})

# ------------------------------------------------------------ post-fit

test_that("simulate() draws correlated residuals", {
  d <- ac_sim(seed = 11, G = 30, K = 5, sd = 1.5)
  f <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(), data = d)
  rho <- autocor_matrix(f)[1, 2]
  s <- as.matrix(simulate(f, nsim = 200, seed = 1))
  r <- array(s - fitted(f), dim = c(5, 30, 200))
  lag1 <- mean(apply(r, 3, function(m) {
    stats::cor(as.vector(m[-5, ]), as.vector(m[-1, ]))
  }))
  expect_equal(lag1, rho, tolerance = 0.05)
  expect_equal(stats::sd(as.vector(r)), sigma(f), tolerance = 0.05)
})

test_that("autocor_matrix() is the fitted correlation and predict is unchanged", {
  d <- ac_sim(seed = 11, G = 30, K = 5)
  f <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(), data = d)
  R <- autocor_matrix(f)
  expect_equal(dim(R), c(5L, 5L))
  expect_equal(diag(R), rep(1, 5), ignore_attr = TRUE)
  expect_equal(R[1, 3], R[1, 2]^2, tolerance = 1e-10)
  expect_identical(dimnames(R)[[1]], as.character(1:5))
  # the mean structure is untouched
  expect_equal(max(abs(fitted(f) - predict(f, type = "response"))), 0)
  nd <- d[1:5, ]
  expect_length(predict(f, newdata = nd), 5L)
  expect_true(all(is.finite(predict(f, newdata = nd,
                                    se.fit = TRUE)$se.fit)))
  # pearson divides by the marginal SD, which is sigma (R is unit
  # diagonal), so it stays the plain standardization
  expect_equal(residuals(f, type = "pearson"),
               residuals(f) / sigma(f), tolerance = 1e-12)
  expect_null(autocor_matrix(frm(bf(y ~ x) + gaussian(), data = d)))
})

test_that("the parameters reach summary, confint and hypothesis", {
  d <- ac_sim(seed = 11, G = 30, K = 5)
  f <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(), data = d)
  expect_true("thetaac_1" %in% rownames(confint(f)))
  vc <- confint_varcorr(f)
  expect_true("ar[1]" %in% vc$term)
  expect_equal(vc$estimate[vc$term == "ar[1]"], autocor_matrix(f)[1, 2],
               tolerance = 1e-10)
  expect_true(all(vc$lwr[vc$term == "ar[1]"] < vc$estimate[vc$term == "ar[1]"]))
  s <- summary(f)
  expect_true("ar[1]" %in% rownames(s$autocor))
  expect_true("ar1" %in% variables(f))
  h <- hypothesis(f, "ar1 - 0.5 = 0")
  expect_true(is.finite(h$se))
  # the fitted R is what the natural-scale name reports
  expect_equal(hypothesis(f, "ar1")$estimate, autocor_matrix(f)[1, 2],
               tolerance = 1e-8)
})

test_that("cosy and unstr name their parameters as brms does", {
  set.seed(5)
  d <- expand.grid(week = 1:4, subj = factor(1:30))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(30, 0, 1.2)[d$subj] +
    stats::rnorm(nrow(d))
  fc <- frm(bf(y ~ x + cosy(week, subj)) + gaussian(), data = d)
  expect_true("cosy" %in% confint_varcorr(fc)$term)
  expect_true("cosy" %in% variables(fc))
  fu <- suppressWarnings(frm(bf(y ~ x + unstr(week, subj)) + gaussian(),
                             data = d))
  tm <- confint_varcorr(fu)$term
  expect_true(all(c("cortime__1__2", "cortime__3__4") %in% tm))
  expect_length(fu$estimates$thetaac, 6L)
  expect_equal(autocor_matrix(fu)[1, 2],
               confint_varcorr(fu)$estimate[tm == "cortime__1__2"],
               tolerance = 1e-10)
})

# --------------------------------------------------------------- structure

test_that("the frame groups rows into time patterns", {
  d <- ac_sim(seed = 2, G = 6, K = 4)
  fr <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(),
            data = d, dry_run = "frame")
  ac <- fr$autocor[[1L]]
  expect_equal(ac$d, 4L)
  expect_equal(ac$n_groups, 6L)
  expect_length(ac$patterns, 1L)          # balanced: one pattern
  expect_equal(ac$patterns[[1L]]$k, 4L)
  expect_equal(ac$patterns[[1L]]$G, 6L)
  expect_equal(ac$npar, 1L)
  expect_length(fr$par_template$thetaac, 1L)
  # every row belongs to exactly one pattern, once
  rows <- sort(unlist(lapply(ac$patterns, `[[`, "rows")))
  expect_equal(rows, seq_len(nrow(d)))
})

test_that("time levels drive the lag; row order is the fallback", {
  d <- ac_sim(seed = 2, G = 6, K = 4)
  # already sorted by (subj, week), so naming the time variable and
  # leaving it out must give the same fit
  f1 <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(), data = d)
  f2 <- frm(bf(y ~ x + ar(gr = subj, cov = TRUE)) + gaussian(), data = d)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-8)
  # shuffling the rows must not change the fit when time is named
  set.seed(4)
  ds <- d[sample(nrow(d)), ]
  f3 <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(), data = ds)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f3)),
               tolerance = 1e-8)
})

test_that("parameter counts follow the structure", {
  d <- ac_sim(seed = 2, G = 8, K = 5)
  np <- function(tm) {
    length(frm(stats::as.formula(paste("y ~ x +", tm)), data = d,
               family = gaussian(),
               dry_run = "frame")$par_template$thetaac)
  }
  expect_equal(np("ar(week, subj, cov = TRUE)"), 1L)
  expect_equal(np("ar(week, subj, p = 3, cov = TRUE)"), 3L)
  expect_equal(np("ma(week, subj, q = 2, cov = TRUE)"), 2L)
  expect_equal(np("arma(week, subj, p = 2, q = 1, cov = TRUE)"), 3L)
  expect_equal(np("cosy(week, subj)"), 1L)
  expect_equal(np("unstr(week, subj)"), 10L)
})

test_that("a non-consecutive whole-number time grid warns", {
  d <- ac_sim(seed = 2, G = 6, K = 4)
  d$week <- c(1, 2, 3, 7)[d$week]
  expect_warning(
    frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(), data = d,
        dry_run = "frame"),
    "not consecutive")
})

# ------------------------------------------------------------------ guards

test_that("only gaussian and student are accepted", {
  d <- ac_sim(seed = 2, G = 8, K = 4)
  d$cnt <- stats::rpois(nrow(d), 3)
  expect_error(
    frm(bf(cnt ~ x + ar(week, subj, cov = TRUE)) + poisson(), data = d),
    "needs a family with real residuals")
  # the refusal names the random-effect spelling of the same idea
  expect_error(
    frm(bf(cnt ~ x + ar(week, subj, cov = TRUE)) + poisson(), data = d),
    "ar1\\(")
  expect_s3_class(
    frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + student(), data = d,
        dry_run = "frame"), "frmtmb_frame")
  expect_error(
    frm(bf(y ~ x + ar(week, subj, cov = TRUE), nu ~ x) + student(),
        data = d),
    "needs a constant nu")
})

test_that("brms's cov = FALSE spelling is refused, not reinterpreted", {
  d <- ac_sim(seed = 2, G = 8, K = 4)
  expect_error(frm(bf(y ~ x + ar(week, subj)) + gaussian(), data = d),
               "needs cov = TRUE")
  expect_error(frm(bf(y ~ x + arma(week, subj)) + gaussian(), data = d),
               "needs cov = TRUE")
  # cosy() and unstr() have no cov argument in brms either
  expect_s3_class(frm(bf(y ~ x + cosy(week, subj)) + gaussian(),
                      data = d, dry_run = "frame"), "frmtmb_frame")
})

test_that("aterms that need a per-row density are refused", {
  d <- ac_sim(seed = 2, G = 8, K = 4)
  d$w <- 1
  d$cc <- 0
  d$sd <- 0.5
  for (lhs in c("y | weights(w)", "y | cens(cc)", "y | trunc(lb = -99)",
                "y | se(sd)")) {
    expect_error(
      frm(bf(stats::as.formula(
        paste(lhs, "~ x + ar(week, subj, cov = TRUE)"))) + gaussian(),
        data = d),
      "cannot be combined with a residual correlation term")
  }
})

test_that("structural combinations are refused", {
  d <- ac_sim(seed = 2, G = 8, K = 4)
  d$y2 <- stats::rnorm(nrow(d))
  expect_error(
    frm(bf(y ~ x, sigma ~ x + ar(week, subj, cov = TRUE)) + gaussian(),
        data = d),
    "can only be written on 'mu'")
  expect_error(
    frm(bf(y ~ x + ar(week, subj, cov = TRUE) + cosy(week, subj)) +
          gaussian(), data = d),
    "residual correlation terms")
  expect_error(
    frm(mvbf(bf(y ~ x + ar(week, subj, cov = TRUE)), bf(y2 ~ x),
             rescor = TRUE) + gaussian(), data = d),
    "rescor = TRUE")
  expect_error(
    frm(bf(y ~ x + ar(week, subj, cov = TRUE)) +
          mixture(gaussian, gaussian), data = d),
    "mixture component")
  expect_error(
    frm(bf(y ~ a * x + ar(week, subj, cov = TRUE), a ~ 1, nl = TRUE) +
          gaussian(), data = d),
    "nonlinear")
  expect_error(
    frm(bf(y ~ x + (1 | subj) + ar(week, subj, cov = TRUE)) + gaussian(),
        data = d, quadrature = TRUE),
    "quadrature = TRUE cannot be combined with the residual")
})

test_that("post-fit paths that need a per-row density are refused", {
  d <- ac_sim(seed = 2, G = 8, K = 4)
  f <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(), data = d)
  expect_error(residuals(f, type = "osa"), "not available")
  # frm_simulate() DOES draw the correlated residual since v0.36: the
  # de novo shim carries the frame's autocor block, and the simulator
  # contract reads it there exactly as simulate() does
  s <- frm_simulate(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(),
                    data = d, nsim = 2, seed = 1,
                    newparams = list(beta = c(0, 1),
                                     betad = as.numeric(f$estimates$betad),
                                     thetaac = as.numeric(f$estimates$thetaac)))
  expect_equal(dim(s), c(nrow(d), 2L))
})

test_that("the time / group argument order is brms's, and says so", {
  d <- ac_sim(seed = 2, G = 8, K = 4)
  # cosy(subj) reads subj as TIME (brms's positional order), so every
  # group has repeated time points
  expect_error(frm(bf(y ~ x + cosy(subj)) + gaussian(), data = d),
               "write cosy\\(gr = ")
  expect_error(frm(bf(y ~ x + unstr(week)) + gaussian(), data = d),
               "needs both a time variable and a grouping variable")
  expect_s3_class(frm(bf(y ~ x + cosy(gr = subj)) + gaussian(), data = d,
                      dry_run = "frame"), "frmtmb_frame")
})

test_that("degenerate and oversized blocks are refused", {
  d <- ac_sim(seed = 2, G = 8, K = 4)
  d$one <- 1
  expect_error(frm(bf(y ~ x + ar(one, subj, cov = TRUE)) + gaussian(),
                   data = d), "at least 2 time points")
  # no grouping variable and no time: the whole data set is one series
  big <- data.frame(t = seq_len(400), y = stats::rnorm(400),
                    x = stats::rnorm(400))
  expect_error(frm(bf(y ~ x + ar(t, cov = TRUE)) + gaussian(),
                   data = big), "cap 300")
  wide <- data.frame(t = seq_len(60), g = factor(1), y = stats::rnorm(60),
                     x = stats::rnorm(60))
  expect_error(frm(bf(y ~ x + unstr(t, g)) + gaussian(), data = wide),
               "cap 50 levels")
})

test_that("REML keeps the correlation parameters in the outer problem", {
  d <- ac_sim(seed = 2, G = 12, K = 4)
  f <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(), data = d,
           REML = TRUE)
  expect_true("thetaac" %in% names(f$opt$par))
  expect_false("beta" %in% names(f$opt$par))
  expect_true("thetaac_1" %in% outer_par_names(f))
  expect_true("thetaac_1" %in% rownames(confint(f)))
})

test_that("an autocor term inside an interaction is refused by name", {
  d <- ac_sim(seed = 3, G = 6, K = 4)
  expect_error(
    frm(bf(y ~ x:ar(week, subj, cov = TRUE)) + gaussian(), data = d,
        dry_run = "frame"),
    "cannot be crossed"
  )
  expect_error(
    frm(bf(y ~ x * cosy(week, subj)) + gaussian(), data = d,
        dry_run = "frame"),
    "separate term"
  )
})
