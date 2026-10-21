# dirichlet_multinomial(): brms 2.23.0's overdispersed multinomial on
# multinomial()'s response layout. The references are an exact form of
# the density (rising factorials as sums of logs, no lgamma difference),
# beta_binomial() at two categories, multinomial() as phi grows, and
# the analytic moments. brms's own program is the gated
# test-nb2dm-brms.R.

dm_data <- function(seed, n = 300, phi = 6, re = FALSE) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                   g = gl(15, n / 15))
  u <- if (re) stats::rnorm(15, 0, 0.5)[dd$g] else 0
  E <- cbind(0, 0.4 + 0.6 * dd$x + u, -0.3 - 0.5 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  dd$N <- sample(5:30, n, TRUE)
  Y <- t(vapply(seq_len(n), function(i) {
    p <- stats::rgamma(3, P[i, ] * phi)
    as.numeric(stats::rmultinom(1, dd$N[i], p / sum(p)))
  }, numeric(3)))
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd
}

# log P(y) with the rising factorials summed term by term: no lgamma
# difference, so it is exact to rounding at any precision
dm_exact <- function(y, alpha) {
  phi <- sum(alpha)
  tot <- sum(y)
  out <- lgamma(tot + 1) - sum(lgamma(y + 1))
  for (k in seq_along(y)) {
    if (y[k] > 0) out <- out + sum(log(alpha[k] + seq_len(y[k]) - 1))
  }
  if (tot > 0) out <- out - sum(log(phi + seq_len(tot) - 1))
  out
}

dm_row <- function(fam, y, eta, phi) {
  dpn <- fam$primary_dpars
  dp <- c(stats::setNames(as.list(eta), dpn), list(phi = phi))
  ym <- matrix(y, 1)
  as.numeric(row_lpdf(fam, ym, ym, dp, list(trials = sum(y)), NULL))
}

test_that("dirichlet_multinomial() is brms's family", {
  fam <- dirichlet_multinomial()
  expect_s3_class(fam, "frmtmb_family")
  expect_identical(fam$family, "dirichlet_multinomial")
  expect_identical(fam$link, "logit")
  expect_identical(fam$link_phi, "log")
  expect_identical(dirichlet_multinomial(link_phi = "softplus")$link_phi,
                   "softplus")
  expect_error(dirichlet_multinomial("probit"),
               paste("'probit' is not a supported link for family",
                     "'dirichlet_multinomial'. Supported links are: 'logit'"),
               fixed = TRUE)
  expect_error(dirichlet_multinomial(link_phi = "logit"),
               "not a supported link for parameter 'phi'", fixed = TRUE)
  f3 <- dirichlet_multinomial(K = 3)
  expect_identical(f3$dpars, c("mu2", "mu3", "phi"))
  expect_identical(f3$primary_dpars, c("mu2", "mu3"))
  expect_identical(brmsfamily("dirichlet_multinomial", K = 3)$dpars,
                   f3$dpars)
})

test_that("a brms dirichlet_multinomial object is read with its refcat", {
  skip_if_not_installed("brms")
  dd <- dm_data(51)
  f <- frm(Y | trials(N) ~ x,
           family = brms::dirichlet_multinomial(refcat = "b"), data = dd)
  expect_identical(names(fixef_by_dpar(f))[1:2], c("mua", "muc"))
})

test_that("the density is the exact Dirichlet-multinomial", {
  set.seed(52)
  worst <- 0
  for (K in c(2L, 3L, 5L)) {
    fam <- dirichlet_multinomial(K = K)
    for (phi in c(1e-3, 0.3, 4, 150, 1e4)) {
      for (r in 1:4) {
        tot <- sample(c(0L, 1L, 9L, 60L), 1)
        eta <- stats::rnorm(K - 1, 0, 1.5)
        p <- exp(c(0, eta)) / sum(exp(c(0, eta)))
        y <- as.vector(stats::rmultinom(1, tot, p))
        a <- dm_row(fam, y, eta, phi)
        b <- dm_exact(y, phi * p)
        # relative to the magnitude of the terms the lgamma form adds,
        # whose rounding bounds what any lgamma evaluation can reach
        scale <- abs(lgamma(phi)) + abs(lgamma(tot + phi)) +
          sum(abs(lgamma(y + phi * p))) + sum(abs(lgamma(phi * p))) +
          abs(lgamma(tot + 1)) + 1
        worst <- max(worst, abs(a - b) / (.Machine$double.eps * scale))
      }
    }
  }
  expect_lt(worst, 64)
})

test_that("the density keeps its digits as phi grows", {
  # phi runs off on data with no overdispersion; brms's lgamma
  # differences lose 2.5e-7 of the density at phi = 1e8 and all of it
  # by 1e16, where the multinomial limit is an ordinary number
  fam <- dirichlet_multinomial(K = 3)
  y <- c(4, 1, 7)
  eta <- c(-0.5, 0.8)
  p <- exp(c(0, eta)) / sum(exp(c(0, eta)))
  for (phi in c(1e6, 1e9, 1e12, 1e16)) {
    a <- dm_row(fam, y, eta, phi)
    b <- dm_exact(y, phi * p)
    # in units of the rounding of the multinomial coefficient, the
    # largest term the density adds at any phi
    scale <- lgamma(sum(y) + 1) + sum(lgamma(y + 1)) + abs(b)
    expect_lt(abs(a - b) / (.Machine$double.eps * scale), 64,
              label = paste("phi", phi))
  }
})

test_that("two categories are beta_binomial(), and a large phi multinomial()",
          {
  bb <- beta_binomial()
  f2 <- dirichlet_multinomial(K = 2)
  for (pt in list(c(y = 3, n = 10, eta = 0.4, phi = 2),
                  c(y = 0, n = 25, eta = -1.2, phi = 40),
                  c(y = 25, n = 25, eta = 2, phi = 0.5))) {
    a <- dm_row(f2, c(pt[["n"]] - pt[["y"]], pt[["y"]]), pt[["eta"]],
                pt[["phi"]])
    b <- as.numeric(row_lpdf(bb, pt[["y"]], pt[["y"]],
                             list(mu = stats::plogis(pt[["eta"]]),
                                  phi = pt[["phi"]]),
                             list(trials = pt[["n"]]), NULL))
    expect_lt(abs(a - b), 1e-10 * abs(b))
  }
  # the multinomial is the limit; the gap shrinks as 1 / phi
  f3 <- dirichlet_multinomial(K = 3)
  mn <- multinomial(K = 3)
  y <- c(4, 1, 7)
  eta <- c(-0.5, 0.8)
  lm <- as.numeric(row_lpdf(mn, matrix(y, 1), matrix(y, 1),
                            list(mu2 = eta[1], mu3 = eta[2]),
                            list(trials = 12), NULL))
  gap <- vapply(c(1e2, 1e3, 1e4), function(ph) {
    abs(dm_row(f3, y, eta, ph) - lm)
  }, 0)
  expect_lt(abs(gap[2] / gap[1] - 0.1), 0.02)
  expect_lt(abs(gap[3] / gap[2] - 0.1), 0.02)
})

test_that("a fit reports phi, brms's fitted() shape and multinomial's mean", {
  dd <- dm_data(53)
  fit <- frm(Y | trials(N) ~ x, family = dirichlet_multinomial(), data = dd)
  mn <- frm(Y | trials(N) ~ x, family = multinomial(), data = dd)
  expect_identical(names(fixef_by_dpar(fit)), c("mub", "muc", "phi"))
  # overdispersed data: the Dirichlet layer is worth many log units
  expect_gt(as.numeric(logLik(fit)) - as.numeric(logLik(mn)),
            0.05 * abs(as.numeric(logLik(mn))))
  ft <- fitted(fit)
  expect_identical(dim(ft), c(nrow(dd), 4L, 3L))
  expect_identical(dimnames(ft)[[3]], c("P(Y = a)", "P(Y = b)", "P(Y = c)"))
  # the expected counts sum to each row's trials
  tot <- apply(ft[, "Estimate", ], 1, sum)
  expect_lt(max(abs(tot / dd$N - 1)), 1e-12)
  # and are the softmax of the predictors times the trials, phi aside
  b <- fixef_by_dpar(fit)
  e2 <- b$mub[[1]] + b$mub[[2]] * dd$x
  e3 <- b$muc[[1]] + b$muc[[2]] * dd$x
  den <- 1 + exp(e2) + exp(e3)
  expect_lt(max(abs(ft[, "Estimate", 2] / (dd$N * exp(e2) / den) - 1)),
            1e-10)
  # a change of reference category changes the coefficients only
  f2 <- frm(Y | trials(N) ~ x, family = dirichlet_multinomial(refcat = "b"),
            data = dd)
  l1 <- as.numeric(logLik(fit))
  expect_lt(abs(as.numeric(logLik(f2)) - l1),
            sqrt(.Machine$double.eps) * abs(l1))
  # a formula on phi, and a group term
  f3 <- frm(bf(Y | trials(N) ~ x, phi ~ z), family = dirichlet_multinomial(),
            data = dd)
  expect_named(fixef_by_dpar(f3)$phi, c("(Intercept)", "z"))
  f4 <- frm(Y | trials(N) ~ x + (1 | g), family = dirichlet_multinomial(),
            data = dm_data(54, re = TRUE))
  expect_identical(f4$opt$convergence, 0L)
})

test_that("on multinomial data phi runs off and the fit stays below the limit",
          {
  # seed 4 of dev/nb2dm-dm-limit.R: with brms's lgamma differences this
  # fit stopped at "false convergence (8)" 1.1e-4 ABOVE the multinomial
  # supremum, which the model can only approach
  set.seed(4)
  n <- 300
  d <- data.frame(x = stats::rnorm(n))
  E <- cbind(0, 0.3 + 0.5 * d$x, -0.2 - 0.4 * d$x)
  P <- exp(E) / rowSums(exp(E))
  d$N <- sample(5:40, n, TRUE)
  Y <- t(vapply(seq_len(n), function(i) {
    as.numeric(stats::rmultinom(1, d$N[i], P[i, ]))
  }, numeric(3)))
  colnames(Y) <- c("a", "b", "c")
  d$Y <- Y
  mn <- frm(Y | trials(N) ~ x, family = multinomial(), data = d)
  dm <- suppressMessages(frm(Y | trials(N) ~ x,
                             family = dirichlet_multinomial(), data = d))
  expect_identical(dm$opt$convergence, 0L)
  lm <- as.numeric(logLik(mn))
  expect_lte(as.numeric(logLik(dm)), lm + .Machine$double.eps * abs(lm))
  # phi has run off (the limit), and the categories' coefficients are
  # the multinomial's
  expect_gt(dm$estimates$betad[[1]], 12)
  b1 <- unlist(fixef_by_dpar(dm)[c("mub", "muc")])
  b0 <- unlist(fixef_by_dpar(mn))
  expect_lt(max(abs(b1 - b0) / abs(b0)), 1e-4)
})

test_that("trials() is required and checked, as for multinomial()", {
  dd <- dm_data(55)
  expect_error(frm(Y ~ x, family = dirichlet_multinomial(), data = dd),
               "Specifying 'trials' is required for this model",
               fixed = TRUE)
  bad <- dd
  bad$N[3] <- bad$N[3] + 1L
  expect_error(frm(Y | trials(N) ~ x, family = dirichlet_multinomial(),
                   data = bad),
               "Number of trials does not match the number of events",
               fixed = TRUE)
  # the refusal names the family, without the bernoulli() hint a
  # binomial gets
  m <- tryCatch(frm(Y ~ x, family = dirichlet_multinomial(), data = dd),
                error = conditionMessage)
  expect_match(m, "dirichlet_multinomial", fixed = TRUE)
  expect_false(grepl("bernoulli", m, fixed = TRUE))
})

test_that("predictions read trials() from newdata", {
  dd <- dm_data(56)
  fit <- frm(Y | trials(N) ~ x, family = dirichlet_multinomial(), data = dd)
  nd <- data.frame(x = c(-1, 0, 1), N = c(4L, 10L, 40L))
  ft <- fitted(fit, newdata = nd)
  expect_lt(max(abs(apply(ft[, "Estimate", ], 1, sum) / nd$N - 1)), 1e-12)
  set.seed(3)
  pr <- predict(fit, newdata = nd, ndraws = 200)
  expect_identical(dim(pr), c(3L, 4L, 3L))
  expect_error(fitted(fit, newdata = data.frame(x = 0)),
               "trials(N) could not be evaluated on newdata", fixed = TRUE)
  s <- simulate(fit, nsim = 2, seed = 5)
  expect_true(all(vapply(s, function(m) all(rowSums(m) == dd$N), NA)))
})

test_that("fitted() holds the Wald bounds of a count between 0 and trials", {
  # the count clamp, not the probability clamp: at x = 8 category b's
  # upper bound would pass the row's 50 trials and is held there, and
  # every bound of b is a count far above 1
  dd <- dm_data(56)
  fit <- frm(Y | trials(N) ~ x, family = dirichlet_multinomial(), data = dd)
  nd <- data.frame(x = c(3, 5, 8), N = c(50L, 50L, 50L))
  ft <- fitted(fit, newdata = nd)
  # held at the trials themselves, not at the estimates' sum, which
  # matches them only to rounding and passed 50 on a Windows runner
  expect_identical(ft[3, "Q97.5", 2], 50)
  expect_true(all(ft[, "Q2.5", 2] > 1))
  expect_true(all(ft[, "Q2.5", ] >= 0) && all(ft[, "Q97.5", ] <= 50))
})

test_that("the simulator's moments are the Dirichlet-multinomial's", {
  fam <- dirichlet_multinomial(K = 3)
  set.seed(57)
  n <- 1e5
  z <- stats::qnorm(1 - 1e-6 / 2)
  for (pt in list(list(eta = c(0.5, -0.4), phi = 2.5, tot = 12),
                  list(eta = c(-1, 1.5), phi = 40, tot = 30))) {
    dp <- list(mu2 = rep(pt$eta[1], n), mu3 = rep(pt$eta[2], n),
               phi = rep(pt$phi, n))
    x <- fam$sim(dp, list(trials = rep(pt$tot, n)), n)
    p <- exp(c(0, pt$eta)) / sum(exp(c(0, pt$eta)))
    infl <- (pt$tot + pt$phi) / (1 + pt$phi)
    for (k in 1:3) {
      m <- pt$tot * p[k]
      v <- pt$tot * p[k] * (1 - p[k]) * infl
      lab <- paste0("phi ", pt$phi, ", category ", k)
      expect_lt(abs(mean(x[, k]) - m), z * sqrt(v / n),
                label = paste(lab, "mean"))
      se_v <- stats::sd((x[, k] - mean(x[, k]))^2) / sqrt(n)
      expect_lt(abs(stats::var(x[, k]) - v), z * se_v,
                label = paste(lab, "variance"))
    }
    # and the covariance of two categories, -T p_j p_k times the same
    # inflation
    cv <- -pt$tot * p[2] * p[3] * infl
    pr <- (x[, 2] - mean(x[, 2])) * (x[, 3] - mean(x[, 3]))
    expect_lt(abs(mean(pr) - cv), z * stats::sd(pr) / sqrt(n),
              label = paste("phi", pt$phi, "covariance"))
  }
})

test_that("the polytomous refusals name dirichlet_multinomial", {
  dd <- dm_data(58)
  fit <- frm(Y | trials(N) ~ x, family = dirichlet_multinomial(), data = dd)
  expect_error(residuals(fit), "'dirichlet_multinomial' family, as in brms",
               fixed = TRUE)
  expect_error(dharma_residuals(fit),
               "no meaning for the 'dirichlet_multinomial' family",
               fixed = TRUE)
  expect_true(frmtmb:::fam_is_polytomous(fit$spec$responses[[1]]$family))
  expect_identical(frmtmb:::mu_stack_dpars(fit$spec$responses[[1]]$family),
                   c("mub", "muc"))
})

test_that("the compatibility table answers for dirichlet_multinomial", {
  tab <- frm_compat("dirichlet_multinomial")
  side <- function(f) tab[tab$feature_a == f | tab$feature_b == f, "status"]
  expect_identical(side("trials()"), "conditional")
  expect_identical(side("simulate"), "works")
  expect_identical(side("fitted"), "works")
  expect_identical(side("cens()"), frm_compat("multinomial")[
    frm_compat("multinomial")$feature_b == "cens()", "status"])
})
