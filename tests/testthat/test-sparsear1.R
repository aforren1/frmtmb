# The O(d) AR(1) density, held against the dense one it replaced.
#
# The dense computation is reproduced INLINE here rather than kept in
# the package: it is the reference, not a code path. It is a verbatim
# copy of what `covstruct_registry$ar1$nll` and
# `covstruct_registry$hetar1$nll` did through v0.42.0, so the gates
# below are the compatibility claim: every fit made before the change
# sees the same number afterwards.

dense_ar1_nll <- function(b, theta, blk) {
  "[<-" <- RTMB::ADoverload("[<-")
  d <- blk$dim
  sd1 <- exp(theta[1])
  rho <- theta[2] / sqrt(1 + theta[2]^2)
  pows <- rep(rho, d)
  pows[1] <- 1
  for (k in seq_len(d - 1L) + 1L) pows[k] <- pows[k - 1L] * rho
  M <- abs(outer(seq_len(d), seq_len(d), "-")) + 1L
  C <- RTMB::matrix(pows[as.vector(M)], d, d)
  Sigma <- sd1^2 * C
  dim(b) <- c(d, length(b) %/% d)
  sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
}

dense_hetar1_nll <- function(b, theta, blk) {
  "[<-" <- RTMB::ADoverload("[<-")
  d <- blk$dim
  sdv <- exp(theta[seq_len(d)])
  rho <- theta[d + 1L] / sqrt(1 + theta[d + 1L]^2)
  pows <- rep(rho, d)
  pows[1] <- 1
  for (k in seq_len(d - 1L) + 1L) pows[k] <- pows[k - 1L] * rho
  M <- abs(outer(seq_len(d), seq_len(d), "-")) + 1L
  C <- RTMB::matrix(pows[as.vector(M)], d, d)
  Sigma <- C * (RTMB::matrix(sdv, ncol = 1) %*% RTMB::matrix(sdv, nrow = 1))
  dim(b) <- c(d, length(b) %/% d)
  sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
}

# 20 configurations of (dim, n_levels, theta, b), drawn once and shared
# by the value gate and the gradient gate.
sparsear1_configs <- function(seed = 20260903, n = 20) {
  set.seed(seed)
  lapply(seq_len(n), function(i) {
    d <- sample(2:24, 1)
    n_levels <- sample(1:12, 1)
    list(dim = d, n_levels = n_levels,
         b = stats::rnorm(d * n_levels),
         # a wide theta: |rho| up to 0.9999 at the tails of this range
         theta_ar1 = c(stats::rnorm(1, 0, 1.2), stats::rnorm(1, 0, 25)),
         theta_het = c(stats::rnorm(d, 0, 1.2), stats::rnorm(1, 0, 25)))
  })
}

test_that("the O(d) ar1 density equals the dense one it replaced", {
  cfgs <- sparsear1_configs()
  worst <- 0
  for (cf in cfgs) {
    blk <- list(dim = cf$dim, n_levels = cf$n_levels)
    got <- covstruct_registry$ar1$nll(cf$b, cf$theta_ar1, blk)
    want <- dense_ar1_nll(cf$b, cf$theta_ar1, blk)
    worst <- max(worst, abs(got - want) / max(1, abs(want)))
  }
  expect_lt(worst, 1e-10)
})

test_that("the O(d) hetar1 density equals the dense one it replaced", {
  cfgs <- sparsear1_configs()
  worst <- 0
  for (cf in cfgs) {
    blk <- list(dim = cf$dim, n_levels = cf$n_levels)
    got <- covstruct_registry$hetar1$nll(cf$b, cf$theta_het, blk)
    want <- dense_hetar1_nll(cf$b, cf$theta_het, blk)
    worst <- max(worst, abs(got - want) / max(1, abs(want)))
  }
  expect_lt(worst, 1e-10)
})

test_that("the gradient through MakeTape agrees with the dense one", {
  cfgs <- sparsear1_configs()
  worst_ar1 <- 0
  worst_het <- 0
  for (cf in cfgs) {
    blk <- list(dim = cf$dim, n_levels = cf$n_levels)
    b <- cf$b
    tp_new <- RTMB::MakeTape(
      function(th) -covstruct_registry$ar1$nll(b, th, blk), cf$theta_ar1)
    tp_old <- RTMB::MakeTape(
      function(th) -dense_ar1_nll(b, th, blk), cf$theta_ar1)
    g1 <- as.vector(tp_new$jacobian(cf$theta_ar1))
    g0 <- as.vector(tp_old$jacobian(cf$theta_ar1))
    worst_ar1 <- max(worst_ar1, max(abs(g1 - g0) / pmax(1, abs(g0))))

    tp_new <- RTMB::MakeTape(
      function(th) -covstruct_registry$hetar1$nll(b, th, blk), cf$theta_het)
    tp_old <- RTMB::MakeTape(
      function(th) -dense_hetar1_nll(b, th, blk), cf$theta_het)
    g1 <- as.vector(tp_new$jacobian(cf$theta_het))
    g0 <- as.vector(tp_old$jacobian(cf$theta_het))
    worst_het <- max(worst_het, max(abs(g1 - g0) / pmax(1, abs(g0))))
  }
  expect_lt(worst_ar1, 1e-8)
  expect_lt(worst_het, 1e-8)
})

test_that("the scaling matches RTMB::dautoreg, which fixes the marginal sd", {
  # An independent read on the derivation in R/covstruct.R: dautoreg's
  # field has unit MARGINAL variance, so `scale =` the marginal standard
  # deviation reproduces our parameterization exactly. If the density
  # were written against the INNOVATION sd instead, this gate is the one
  # that would catch it.
  cfgs <- sparsear1_configs(seed = 4242, n = 8)
  worst <- 0
  for (cf in cfgs) {
    blk <- list(dim = cf$dim, n_levels = cf$n_levels)
    m <- matrix(cf$b, cf$dim, cf$n_levels)
    rho <- cf$theta_ar1[2] / sqrt(1 + cf$theta_ar1[2]^2)
    ref <- sum(vapply(seq_len(cf$n_levels), function(j) {
      RTMB::dautoreg(m[, j], phi = rho, log = TRUE,
                     scale = exp(cf$theta_ar1[1]))
    }, numeric(1)))
    got <- covstruct_registry$ar1$nll(cf$b, cf$theta_ar1, blk)
    worst <- max(worst, abs(got - ref) / max(1, abs(ref)))

    rho <- cf$theta_het[cf$dim + 1L] /
      sqrt(1 + cf$theta_het[cf$dim + 1L]^2)
    ref <- sum(vapply(seq_len(cf$n_levels), function(j) {
      RTMB::dautoreg(m[, j], phi = rho, log = TRUE,
                     scale = exp(cf$theta_het[seq_len(cf$dim)]))
    }, numeric(1)))
    got <- covstruct_registry$hetar1$nll(cf$b, cf$theta_het, blk)
    worst <- max(worst, abs(got - ref) / max(1, abs(ref)))
  }
  expect_lt(worst, 1e-10)
})

test_that("the density is the one vcov() and chol_L() describe", {
  # The registry's other accessors were not touched, so they are the
  # third statement of the same covariance: Sigma from vcov(), L L' from
  # chol_L(), and the density itself must all be the same matrix.
  blk <- list(dim = 6L, n_levels = 1L, cnms = paste0("t", 1:6))
  th_ar1 <- c(0.3, -1.4)
  th_het <- c(seq(-0.5, 0.7, length.out = 6), 0.8)
  b <- stats::rnorm(6)
  for (nm in c("ar1", "hetar1")) {
    th <- if (nm == "ar1") th_ar1 else th_het
    V <- covstruct_registry[[nm]]$vcov(th, blk)
    L <- covstruct_registry[[nm]]$chol_L(th, blk)
    expect_equal(unname(as.matrix(L %*% t(L))), unname(V), tolerance = 1e-12)
    expect_equal(covstruct_registry[[nm]]$nll(b, th, blk),
                 as.numeric(RTMB::dmvnorm(b, 0, unname(V), log = TRUE)),
                 tolerance = 1e-10)
  }
})

test_that("the registry contracts around ar1 are unchanged", {
  expect_equal(covstruct_registry$ar1$npar(9L), 2L)
  expect_equal(covstruct_registry$ar1$sd_idx(9L), 1L)
  expect_equal(covstruct_registry$ar1$start(9L), c(0, 0))
  expect_equal(covstruct_registry$ar1$cor_spec(9L),
               list(kind = "ar1", d = 9L, idx = 2L))
  expect_null(covstruct_registry$ar1$cor_spec(1L))
  expect_equal(covstruct_registry$hetar1$npar(9L), 10L)
  expect_equal(covstruct_registry$hetar1$sd_idx(9L), seq_len(9L))
  expect_equal(covstruct_registry$hetar1$start(9L), numeric(10L))
  expect_equal(covstruct_registry$hetar1$cor_spec(9L),
               list(kind = "ar1", d = 9L, idx = 10L))
})

test_that("a long ar1 series now tapes and fits", {
  # d = 400 is past the point where the dense route was usable: the
  # 400 x 400 Cholesky went on the tape and one gradient took seconds.
  skip_on_cran()
  set.seed(11)
  d <- 400L
  n_g <- 3L
  rho <- 0.7
  sd_u <- 0.9
  u <- do.call(c, lapply(seq_len(n_g), function(i) {
    as.vector(stats::arima.sim(list(ar = rho), n = d, sd = sd_u *
                                 sqrt(1 - rho^2)))
  }))
  dd <- data.frame(
    tim = factor(rep(seq_len(d), n_g), levels = seq_len(d)),
    g = factor(rep(seq_len(n_g), each = d)),
    y = 1 + u + stats::rnorm(d * n_g, 0, 0.5)
  )
  t0 <- proc.time()[["elapsed"]]
  fit <- frm(bf(y ~ 1 + ar1(tim + 0 | g)) + gaussian(), data = dd)
  el <- proc.time()[["elapsed"]] - t0
  expect_s3_class(fit, "frmtmb_fit")
  expect_true(is.finite(logLik(fit)))
  # a generous ceiling: the point is that it is seconds, not minutes
  expect_lt(el, 120)
  V <- VarCorr(fit)[[1]]
  expect_lt(abs(V[1, 2] / V[1, 1] - rho), 0.15)
  expect_lt(abs(sqrt(V[1, 1]) - sd_u), 0.4)
})

test_that("sdv_multi's grammar mapping is exact (the acceptance proof)", {
  # The multivariate stochastic volatility model of Skaug and Yu (2014),
  # in miniature. The full n = 945 run lives in dev/tmb-examples-check.R;
  # this is its regression form, at a length that keeps the file quick.
  #
  # The claim is about the two LIKELIHOODS, not about two optimizers, so
  # the reference objective is evaluated at frmtmb's estimates rather
  # than optimized on its own. The reference is written here in the
  # upstream compact spelling.
  skip_on_cran()
  set.seed(7)
  n <- 150L
  p <- 3L
  phi <- c(0.9, 0.8, 0.85)
  sig <- c(0.35, 0.4, 0.3)
  mu_x <- c(-0.5, -0.8, -0.2)
  h <- sapply(seq_len(p), function(j) {
    as.vector(stats::arima.sim(list(ar = phi[j]), n = n, sd = sig[j]))
  })
  R <- matrix(c(1, 0.4, 0.2, 0.4, 1, -0.3, 0.2, -0.3, 1), p, p)
  Lr <- t(chol(R))
  X <- t(vapply(seq_len(n), function(i) {
    exp(0.5 * (mu_x + h[i, ])) * as.vector(Lr %*% stats::rnorm(p))
  }, numeric(p)))

  dd <- data.frame(x1 = X[, 1], x2 = X[, 2], x3 = X[, 3],
                   tim = factor(seq_len(n), levels = seq_len(n)),
                   g = factor(rep(1L, n)))
  fit <- frm(mvbf(bf(x1 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
                  bf(x2 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
                  bf(x3 ~ 0, sigma ~ 1 + ar1(tim + 0 | g)) + gaussian(),
                  rescor = TRUE), data = dd)

  # the reference, in the upstream compact spelling
  ref_f <- function(parms) {
    "[<-" <- RTMB::ADoverload("[<-")
    s <- exp(parms$log_sigma)
    ph <- parms$phi
    sigma_init <- s / sqrt(1 - ph^2)
    hh <- parms$h
    nll <- 0
    for (j in seq_len(p)) {
      nll <- nll - RTMB::dautoreg(hh[, j], phi = ph[j],
                                  scale = sigma_init[j], log = TRUE)
    }
    L <- diag(p)
    L[lower.tri(L)] <- parms$off_diag_x
    # rowSums() strips the advector class; the matrix product does not
    L <- L / as.vector(sqrt((L * L) %*% rep(1, p)))
    Rm <- L %*% t(L)
    mux <- RTMB::matrix(rep(parms$mu_x, each = n), n, p)
    sy <- exp(0.5 * (hh + mux))
    nll - sum(RTMB::dmvnorm(X, 0, Rm, scale = sy, log = TRUE))
  }
  obj <- RTMB::MakeADFun(ref_f, list(phi = rep(0.9, p),
                                     log_sigma = rep(-1, p),
                                     mu_x = rep(-0.5, p),
                                     off_diag_x = rep(0, p),
                                     h = matrix(0, n, p)),
                         random = "h", silent = TRUE)

  # the map: frmtmb's sigma has a log link, so its ar1() block is h / 2
  # and its sigma intercept is mu_x / 2
  th <- fit$estimates$theta
  t2 <- th[seq_len(p) * 2L]
  ph <- t2 / sqrt(1 + t2^2)
  sigma_init <- 2 * exp(th[seq_len(p) * 2L - 1L])
  par_at <- c(ph, log(sigma_init * sqrt(1 - ph^2)),
              2 * fit$estimates$betad, fit$estimates$thetar)
  cross <- obj$fn(par_at)
  expect_lt(abs(cross + as.numeric(logLik(fit))) /
              abs(as.numeric(logLik(fit))), 1e-6)
})
