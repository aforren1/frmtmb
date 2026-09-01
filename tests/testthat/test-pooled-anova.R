# Pooled model comparison over imputations: the D1, D2 and D3 rules on
# frmtmb_multiple. Two independent references are used. Against mice the
# comparison is end to end, and it is exact only where our likelihood and
# covariance coincide with the reference fit's: a poisson GLM has no
# dispersion parameter, so ML beta-hat, the inverse-information vcov and
# df.residual all match glm() exactly. D3 is checked against the
# Meng-Rubin formula written out by hand, because mice::D3() re-estimates
# the nuisance parameters instead of plugging the pooled ones in.

# Deterministic MAR imputations of a single predictor: draw the missing
# values from the observed marginal. Good enough to give the pooling
# rules a nonzero between-imputation variance.
# Relative agreement across the whole row; df2 runs into the thousands,
# where an absolute tolerance says nothing.
expect_rel_equal <- function(a, ref, tol) {
  got <- unlist(a[c("statistic", "df1", "df2", "p", "riv")])
  testthat::expect_lt(max(abs(unname(got) / unname(ref) - 1)), tol)
}

make_imps <- function(d, var, m, seed) {
  set.seed(seed)
  obs <- d[[var]][!is.na(d[[var]])]
  lapply(seq_len(m), function(i) {
    di <- d
    na <- is.na(di[[var]])
    di[[var]][na] <- sample(obs, sum(na), replace = TRUE)
    di
  })
}

test_that("D1 and D2 match mice on a poisson GLM", {
  skip_if_not_installed("mice")
  set.seed(9)
  n <- 200
  dp <- data.frame(x = rnorm(n), z = rnorm(n))
  dp$y <- rpois(n, exp(0.3 + 0.4 * dp$x + 0.2 * dp$z))
  dp$x[runif(n) < 0.2] <- NA
  imp <- mice::mice(dp, m = 8, seed = 3, printFlag = FALSE)
  cmp <- mice::complete(imp, "all")
  # nlminb stops ~4e-06 from the IRLS optimum here, which is the whole
  # residual disagreement below; BFGS to reltol 1e-14 removes it and
  # leaves the pooling algebra as the only thing under test
  ctl <- frmtmb_control(optimizer = "optim",
                        optCtrl = list(method = "BFGS", maxit = 1000,
                                       reltol = 1e-14))
  m1 <- frm_multiple(bf(y ~ x + z) + poisson(), data = cmp, control = ctl)
  m0 <- frm_multiple(bf(y ~ x) + poisson(), data = cmp, control = ctl)
  g1 <- with(imp, glm(y ~ x + z, family = poisson))
  g0 <- with(imp, glm(y ~ x, family = poisson))

  # no dispersion parameter, so mice's dfcom convention is ours too
  expect_equal(df.residual(m1$fits[[1]]),
               df.residual(g1$analyses[[1]]))

  expect_rel_equal(anova(m1, m0, method = "D1"),
                   as.numeric(mice::D1(g1, g0)$result), tol = 1e-5)
  expect_rel_equal(anova(m1, m0, method = "D2", use = "wald"),
                   as.numeric(mice::D2(g1, g0)$result), tol = 1e-5)
  # the likelihood form only needs logLik, which matches glm() to
  # optimizer noise, so it is tighter still
  expect_rel_equal(anova(m1, m0, method = "D2", use = "likelihood"),
                   as.numeric(mice::D2(g1, g0, use = "likelihood")$result),
                   tol = 1e-8)
  # mice's fit0 = NULL default zeroes the slopes and keeps a free
  # intercept, so the constraint set is the slopes
  expect_rel_equal(anova(m1, constraint = c("x", "z"), method = "D1"),
                   as.numeric(mice::D1(g1)$result), tol = 1e-5)
})

test_that("D3 reproduces the Meng-Rubin formula written out by hand", {
  set.seed(73)
  n <- 200
  x <- rnorm(n)
  z <- rnorm(n, 0.4 * x, 1)
  dd <- data.frame(y = rnorm(n, 1 + 0.5 * x + 0.3 * z, 1), x = x, z = z)
  dd$x[runif(n) < 0.25] <- NA
  cmp <- make_imps(dd, "x", 6, seed = 5)
  m1 <- frm_multiple(bf(y ~ x + z) + gaussian(), data = cmp)
  m0 <- frm_multiple(bf(y ~ x) + gaussian(), data = cmp)
  a3 <- anova(m1, m0, method = "D3")

  # reference: pool the optimizer's parameter vector (beta and log
  # sigma), then evaluate each imputation's gaussian log-likelihood
  # there by hand
  dev_L <- function(mf, form) {
    pb <- rowMeans(vapply(mf$fits, function(f) f$opt$par,
                          numeric(length(mf$fits[[1]]$opt$par))))
    nms <- names(mf$fits[[1]]$opt$par)
    b <- pb[nms == "beta"]
    lsig <- pb[nms == "betad"]
    vapply(seq_along(mf$fits), function(i) {
      X <- stats::model.matrix(form, cmp[[i]])
      -2 * sum(stats::dnorm(cmp[[i]]$y, as.vector(X %*% b),
                            exp(lsig), log = TRUE))
    }, 0)
  }
  m <- 6
  k <- 1
  devM <- mean(vapply(m0$fits, function(f) 2 * f$opt$objective, 0) -
                 vapply(m1$fits, function(f) 2 * f$opt$objective, 0))
  devL <- mean(dev_L(m0, ~ x) - dev_L(m1, ~ x + z))
  r <- (m + 1) / (k * (m - 1)) * (devM - devL)
  stat <- devL / (k * (1 + r))
  tdf <- k * (m - 1)
  w <- 4 + (tdf - 4) * (1 + (1 - 2 / tdf) / r)^2

  expect_equal(a3$statistic, stat, tolerance = 1e-8)
  expect_equal(a3$riv, r, tolerance = 1e-8)
  expect_equal(a3$df2, w, tolerance = 1e-8)
  expect_equal(a3$df1, k)
  expect_equal(a3$p, stats::pf(stat, k, w, lower.tail = FALSE),
               tolerance = 1e-8)
  expect_true(a3$riv > 0)

  # evaluating at the pooled vector must not disturb the fits: the
  # objective is restored to its optimum, so the post-fit surface still
  # matches an untouched fit of the same imputation
  f_ind <- frm(bf(y ~ x + z) + gaussian(), data = cmp[[2]])
  expect_vector_equal(diag(vcov(m1$fits[[2]])), diag(vcov(f_ind)),
                      tol = 1e-8)
  expect_equal(as.numeric(logLik(m1$fits[[2]])),
               as.numeric(logLik(f_ind)), tolerance = 1e-8)
})

test_that("identical imputations collapse D3 and D2 to the single LRT", {
  set.seed(41)
  n <- 120
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  d$y <- rnorm(n, 0.5 + 0.4 * d$x + 0.35 * d$z)
  m1 <- frm_multiple(bf(y ~ x + z) + gaussian(), data = list(d, d, d))
  m0 <- frm_multiple(bf(y ~ x) + gaussian(), data = list(d, d, d))
  av <- anova(m1$fits[[1]], m0$fits[[1]])
  lrt <- av$Chisq[2]

  # zero between-imputation variance means zero relative increase in
  # variance, an infinite denominator df and the complete-data test
  for (a in list(anova(m1, m0, method = "D3"),
                 anova(m1, m0, method = "D2", use = "likelihood"))) {
    expect_lt(a$riv, 1e-8)
    expect_gt(a$df2, 1e6)
    expect_equal(a$statistic, lrt / a$df1, tolerance = 1e-8)
    expect_equal(a$p, av$`Pr(>Chisq)`[2], tolerance = 1e-6)
  }
  # the Wald rules collapse to the complete-data Wald test instead
  aw <- anova(m1, m0, method = "D1")
  cf <- m1$fits[[1]]$estimates$beta
  expect_equal(aw$statistic,
               (cf[["z"]] / sqrt(vcov(m1$fits[[1]])["z", "z"]))^2,
               tolerance = 1e-8)
  expect_lt(aw$riv, 1e-8)
})

test_that("D-statistics agree with the per-imputation LRTs on a GLMM", {
  set.seed(62)
  n_g <- 25
  n_per <- 8
  n <- n_g * n_per
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- rnorm(n)
  z <- rnorm(n)
  b0 <- rnorm(n_g, 0, 0.7)
  d <- data.frame(y = 1 + 0.8 * z + b0[g] + rnorm(n, 0, 0.6),
                  x = x, z = z, g = g)
  d$x[sample(n, 40)] <- NA
  cmp <- make_imps(d, "x", 4, seed = 7)
  form1 <- bf(y ~ x + z + (1 | g)) + gaussian()
  form0 <- bf(y ~ x + (1 | g)) + gaussian()
  m1 <- frm_multiple(form1, data = cmp)
  m0 <- frm_multiple(form0, data = cmp)

  # z has a real effect, so every imputation rejects on its own; the
  # pooled rules must agree (no mice reference exists for a GLMM)
  per_p <- vapply(seq_len(4), function(i) {
    anova(m1$fits[[i]], m0$fits[[i]])$`Pr(>Chisq)`[2]
  }, 0)
  expect_true(all(per_p < 0.01))
  for (meth in c("D1", "D2", "D3")) {
    a <- anova(m1, m0, method = meth)
    expect_s3_class(a, "frmtmb_pooled_anova")
    expect_equal(nrow(a), 1L)
    expect_equal(a$df1, 1)
    expect_true(is.finite(a$statistic) && a$statistic > 0)
    expect_true(a$riv >= 0 && is.finite(a$riv))
    expect_true(a$p >= 0 && a$p < 0.01)
  }

  # the same comparison for a predictor with no effect must not reject
  mn1 <- frm_multiple(bf(y ~ x + z + (1 | g)) + gaussian(), data = cmp)
  mn0 <- frm_multiple(bf(y ~ z + (1 | g)) + gaussian(), data = cmp)
  per_p0 <- vapply(seq_len(4), function(i) {
    anova(mn1$fits[[i]], mn0$fits[[i]])$`Pr(>Chisq)`[2]
  }, 0)
  expect_true(all(per_p0 > 0.2))
  for (meth in c("D1", "D2", "D3")) {
    expect_true(anova(mn1, mn0, method = meth)$p > 0.2)
  }

  # D3 also pools the variance components, so a random-effect
  # comparison is available where the Wald rules have nothing to test
  mr0 <- frm_multiple(bf(y ~ x + z) + gaussian(), data = cmp)
  a3 <- anova(m1, mr0, method = "D3")
  expect_equal(a3$df1, 1)
  expect_true(a3$p < 0.01)
  expect_error(anova(m1, mr0, method = "D1"), "covariance parameters")
})

test_that("pooled anova refuses incomparable inputs", {
  set.seed(15)
  n <- 80
  d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n), q = rnorm(n))
  d$y <- rnorm(n, 0.3 + 0.5 * d$x)
  imps <- list(d, d, d)
  m1 <- frm_multiple(bf(y ~ x + z) + gaussian(), data = imps)
  m0 <- frm_multiple(bf(y ~ x) + gaussian(), data = imps)

  expect_error(anova(m1), "second frmtmb_multiple")
  expect_error(anova(m1, m0, m0), "two frmtmb_multiple fits")
  expect_error(anova(m1, constraint = "z", method = "D3"),
               "needs a Wald rule")
  expect_error(anova(m1, constraint = "w", method = "D1"),
               "names no coefficient")
  expect_error(anova(m1, m1, method = "D3"), "different size")

  # different numbers of imputations
  m0b <- frm_multiple(bf(y ~ x) + gaussian(), data = imps[1:2])
  expect_error(anova(m1, m0b), "same imputations")

  # same m, but one imputation lost rows to NA
  dna <- d
  dna$z[1:5] <- NA
  m0c <- frm_multiple(bf(y ~ x + z) + gaussian(),
                      data = list(d, dna, d))
  expect_error(anova(m0c, m0), "same imputed datasets")

  # non-nested coefficient sets (q is absent from the larger model)
  mbig <- frm_multiple(bf(y ~ x + z + w) + gaussian(), data = imps)
  m0d <- frm_multiple(bf(y ~ x + q) + gaussian(), data = imps)
  expect_error(anova(mbig, m0d, method = "D1"), "not nested")

  # REML, as anova.frmtmb_fit refuses it
  mr <- frm_multiple(bf(y ~ x + z) + gaussian(), data = imps, REML = TRUE)
  expect_error(anova(mr, m0), "REML = FALSE")
})

test_that("the pooled anova table prints its models and rule", {
  set.seed(28)
  n <- 60
  d <- data.frame(x = rnorm(n))
  d$y <- rnorm(n, 0.2 + 0.6 * d$x)
  imps <- list(d, d, d)
  m1 <- frm_multiple(bf(y ~ x) + gaussian(), data = imps)
  m0 <- frm_multiple(bf(y ~ 1) + gaussian(), data = imps)
  out <- capture.output(print(anova(m1, m0, method = "D2", use = "wald")))
  expect_true(any(grepl("3 imputations \\(D2, wald\\)", out)))
  expect_true(any(grepl("Model 1: y ~ x", out)))
  expect_true(any(grepl("Model 2: y ~ 1", out)))
  expect_true(any(grepl("statistic", out)))
  out2 <- capture.output(print(anova(m1, constraint = "x",
                                     method = "D1")))
  expect_true(any(grepl("constraint: x = 0", out2)))
})
