# v0.18: mo()/mi() interactions, gp() terms (exact + Hilbert-space),
# group-level latent-class mixtures.

test_that("mo() interactions share the simplex and match direct ML", {
  set.seed(11)
  n <- 500
  inc <- sample(0:3, n, replace = TRUE)
  z <- rnorm(n)
  cz_t <- c(0, 0.5, 0.8, 1)
  y <- 1 + (1.5 + 0.8 * z) * cz_t[inc + 1] + 0.3 * z + rnorm(n, 0, 0.7)
  dd <- data.frame(y = y, inc = inc, z = z)
  fit <- frm(bf(y ~ mo(inc) * z) + gaussian(), data = dd)

  nll <- function(p) {
    zr <- exp(c(0, p[5:7]))
    cz0 <- c(0, cumsum(zr / sum(zr)))
    m <- cz0[inc + 1] * 3
    -sum(stats::dnorm(y, p[1] + p[2] * z + p[3] * m + p[4] * m * z,
                      exp(p[8]), log = TRUE))
  }
  op <- stats::optim(c(1, 0.3, 0.5, 0.3, 0, 0, 0, log(0.7)), nll,
                     method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-6)
  expect_setequal(names(fixef(fit)$mu),
                  c("(Intercept)", "z", "moinc", "moinc:z"))
  # one shared simplex for both mo terms
  expect_length(grep("^zeta", names(fit$frame$par_template)), 1L)
  expect_true(is.finite(predict(fit, newdata = data.frame(inc = 2,
                                                          z = 1))))
})

test_that("mi() interactions stay linear-gaussian exact", {
  set.seed(13)
  n <- 250
  z <- rnorm(n)
  x <- rnorm(n, 0.5 + 0.8 * z, 0.7)
  y <- rnorm(n, 1 + (0.6 + 0.4 * z) * x + 0.3 * z, 0.9)
  x_mis <- x
  mis <- sort(sample(n, 70))
  x_mis[mis] <- NA
  dd <- data.frame(y = y, x = x_mis, z = z)
  fit <- frm(bf(y ~ mi(x) * z) + gaussian() +
               bf(x | mi() ~ z) + gaussian(), data = dd)

  obs <- !is.na(x_mis)
  nll <- function(p) {
    sy <- exp(p[7]); sx <- exp(p[8])
    mux <- p[5] + p[6] * z
    cz <- p[2] + p[4] * z
    -(sum(stats::dnorm(y[obs], p[1] + cz[obs] * x_mis[obs] +
                         p[3] * z[obs], sy, log = TRUE)) +
        sum(stats::dnorm(x_mis[obs], mux[obs], sx, log = TRUE)) +
        sum(stats::dnorm(y[!obs], p[1] + cz[!obs] * mux[!obs] +
                           p[3] * z[!obs],
                         sqrt(sy^2 + cz[!obs]^2 * sx^2), log = TRUE)))
  }
  op <- stats::optim(c(1, 0.6, 0.3, 0.4, 0.5, 0.8, log(0.9), log(0.7)),
                     nll, method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-6)
  expect_true("mix:z" %in% names(fixef(fit)$y_mu))
})

test_that("exact gp() matches direct GP marginal ML", {
  set.seed(17)
  n <- 120
  xg <- round(runif(n, 0, 10), 1)
  pos <- sort(unique(xg))
  Dm <- abs(outer(pos, pos, "-"))
  K <- 1.5^2 * exp(-Dm^2 / 8)
  u <- drop(crossprod(chol(K + diag(1e-8, length(pos))),
                      rnorm(length(pos))))
  dg <- data.frame(y = 2 + u[match(xg, pos)] + rnorm(n, 0, 0.5), x = xg)
  fg <- frm(bf(y ~ gp(x)) + gaussian(), data = dg)

  nll <- function(p) {
    Kp <- exp(2 * p[2]) * exp(-Dm^2 / (2 * exp(2 * p[3])))
    Zi <- outer(xg, pos, `==`) * 1
    S <- Zi %*% Kp %*% t(Zi) + exp(2 * p[4]) * diag(n)
    r <- dg$y - p[1]
    0.5 * (determinant(S)$modulus + sum(r * solve(S, r)) +
             n * log(2 * pi))
  }
  op <- stats::optim(c(2, log(1.5), log(2), log(0.5)), nll,
                     method = "BFGS",
                     control = list(reltol = 1e-12, maxit = 2000))
  expect_lt(abs(as.numeric(logLik(fg)) + op$value), 1e-4)

  # observed-position prediction works; unseen positions krige
  # (closed-form validation in test-gp-multidim.R)
  expect_equal(unname(predict(fg, newdata = dg[1:5, ])),
               unname(fitted(fg)[1:5]), tolerance = 1e-10)
  p_new <- predict(fg, newdata = data.frame(x = 0.05), se.fit = TRUE)
  expect_true(is.finite(p_new$fit) && is.finite(p_new$se.fit))
  cv <- confint_varcorr(fg)
  expect_setequal(cv$term, c("sd(gp)", "range(gp)"))

  # the Hilbert-space approximation converges to the exact answer and
  # predicts anywhere. The boundary convention is brms's: inputs are
  # rescaled by the largest pairwise distance, so L is exactly c (1.25 by
  # default) and k = 40 lands within 5e-3 logLik of the exact fit.
  fh <- frm(bf(y ~ gp(x, k = 40)) + gaussian(), data = dg)
  expect_equal(fh$frame$linpreds[["y.mu"]]$gps[[1]]$L, 1.25,
               tolerance = 1e-12)
  expect_lt(abs(as.numeric(logLik(fh)) - as.numeric(logLik(fg))), 5e-3)
  nd <- data.frame(x = seq(0.5, 9.5, by = 0.5))
  ph <- predict(fh, newdata = nd, se.fit = TRUE)
  expect_true(all(is.finite(ph$fit)) && all(is.finite(ph$se.fit)))
  # curves agree where both are defined
  p_e <- predict(fg, newdata = data.frame(x = pos))
  p_h <- predict(fh, newdata = data.frame(x = pos))
  expect_lt(max(abs(p_e - p_h)), 0.01)

  # newdata rebuilds the basis from the stored scaling, so in-sample rows
  # reproduce the fitted values exactly
  expect_equal(unname(predict(fh, newdata = dg)), unname(fitted(fh)),
               tolerance = 1e-12)

  # the lengthscale is estimated on the rescaled inputs but reported in
  # data units: it agrees with the exact fit's range, and rescaling the
  # covariate rescales the reported range by the same factor
  cvh <- confint_varcorr(fh)
  expect_setequal(cvh$term, c("sd(gp)", "range(gp)"))
  expect_lt(abs(cvh$estimate[cvh$term == "range(gp)"] -
                  cv$estimate[cv$term == "range(gp)"]), 0.01)
  dk <- transform(dg, x = x * 1000)
  fk <- frm(bf(y ~ gp(x, k = 40)) + gaussian(), data = dk)
  expect_equal(as.numeric(logLik(fk)), as.numeric(logLik(fh)),
               tolerance = 1e-8)
  cvk <- confint_varcorr(fk)
  expect_equal(cvk$estimate[cvk$term == "range(gp)"],
               cvh$estimate[cvh$term == "range(gp)"] * 1000,
               tolerance = 1e-6)
})

test_that("group-level latent-class mixtures match direct ML", {
  set.seed(19)
  ng <- 40
  m <- 6
  cls <- rbinom(ng, 1, 0.4)
  g <- rep(seq_len(ng), each = m)
  y <- rnorm(ng * m, ifelse(cls == 1, 2, -1)[g], 0.8)
  dd <- data.frame(y = y, g = factor(g))
  fit <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g),
             data = dd)

  nll <- function(p) {
    pi1 <- stats::plogis(p[5])
    M <- matrix(y, nrow = m)
    l1 <- colSums(stats::dnorm(M, p[1], exp(p[3]), log = TRUE))
    l2 <- colSums(stats::dnorm(M, p[2], exp(p[4]), log = TRUE))
    -sum(log(pi1 * exp(l1) + (1 - pi1) * exp(l2)))
  }
  op <- stats::optim(c(2, -1, log(0.8), log(0.8), 0), nll,
                     method = "BFGS",
                     control = list(reltol = 1e-13, maxit = 5000))
  expect_lt(abs(as.numeric(logLik(fit)) + op$value), 1e-6)

  # posterior class probabilities recover the classes (up to label
  # swap, which mixture likelihoods cannot identify)
  P <- mixture_probs(fit)
  expect_equal(dim(P), c(ng, 2L))
  expect_equal(unname(rowSums(P)), rep(1, ng), tolerance = 1e-10)
  acc <- mean((P[, 1] > 0.5) == (cls == 1))
  expect_gt(max(acc, 1 - acc), 0.95)

  # observation-level mixture_probs work too
  f0 <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = dd)
  P0 <- mixture_probs(f0)
  expect_equal(dim(P0), c(ng * m, 2L))

  # since v0.19 random effects and simulate() work with groups=
  fre <- frm(bf(y ~ 1 + (1 | g)) +
               mixture(gaussian(), gaussian(), groups = ~g), data = dd)
  expect_gte(as.numeric(logLik(fre)), as.numeric(logLik(fit)) - 1e-6)
  expect_equal(nrow(simulate(fit, nsim = 1, seed = 1)), nrow(dd))
})
