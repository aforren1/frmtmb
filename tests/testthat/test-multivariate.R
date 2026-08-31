sim_mv_data <- function(seed = 61, n = 300, n_g = 15, re_cor = 0,
                        res_cor = 0) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_g), length.out = n))
  Su <- matrix(c(0.7^2, re_cor * 0.7 * 0.4,
                 re_cor * 0.7 * 0.4, 0.4^2), 2)
  U <- matrix(rnorm(n_g * 2), n_g) %*% chol(Su)
  Se <- matrix(c(0.8^2, res_cor * 0.8 * 0.5,
                 res_cor * 0.8 * 0.5, 0.5^2), 2)
  E <- matrix(rnorm(n * 2), n) %*% chol(Se)
  x <- rnorm(n)
  data.frame(
    x = x, g = g,
    y1 = 1 + 0.5 * x + U[g, 1] + E[, 1],
    y2 = -1 + 0.2 * x + U[g, 2] + E[, 2]
  )
}

test_that("independent multivariate fit equals the sum of separate fits", {
  dd <- sim_mv_data()
  mv <- frm(mvbf(bf(y1 ~ x + (1 | g)) + gaussian(),
                 bf(y2 ~ x + (1 | g)) + gaussian()), data = dd)
  f1 <- frm(bf(y1 ~ x + (1 | g)) + gaussian(), dd)
  f2 <- frm(bf(y2 ~ x + (1 | g)) + gaussian(), dd)
  expect_lt(abs(as.numeric(logLik(mv)) -
                  (as.numeric(logLik(f1)) + as.numeric(logLik(f2)))),
            1e-6)
  expect_vector_equal(fixef(mv)$y1_mu, fixef(f1)$mu, tol = 1e-5)
  expect_vector_equal(fixef(mv)$y2_mu, fixef(f2)$mu, tol = 1e-5)
})

test_that("bf() + bf() and mvbind() build multivariate formulas", {
  m1 <- bf(y1 ~ x) + gaussian() + bf(y2 ~ x) + gaussian()
  expect_s3_class(m1, "frmtmb_mvformula")
  expect_length(m1$forms, 2)
  m2 <- bf(mvbind(y1, y2) ~ x, family = gaussian())
  expect_s3_class(m2, "frmtmb_mvformula")
  expect_identical(deparse1(m2$forms[[2]]$formula), "y2 ~ x")
})

test_that("rescor matches a hand-rolled RTMB reference", {
  dd <- sim_mv_data(seed = 62, res_cor = 0.5)
  fit <- frm(mvbf(bf(y1 ~ x) + gaussian(),
                  bf(y2 ~ x) + gaussian(), rescor = TRUE), data = dd)

  y1 <- dd$y1; y2 <- dd$y2; x <- dd$x
  nll_ref <- function(p) {
    "c" <- RTMB::ADoverload("c")
    mu1 <- p$b1[1] + p$b1[2] * x
    mu2 <- p$b2[1] + p$b2[2] * x
    s1 <- exp(p$ls1); s2 <- exp(p$ls2)
    rho <- p$tr / sqrt(1 + p$tr^2)
    C <- RTMB::matrix(c(1, rho, rho, 1), 2, 2)
    Z <- RTMB::matrix(c((y1 - mu1) / s1, (y2 - mu2) / s2), length(x), 2)
    -sum(RTMB::dmvnorm(Z, 0, C, log = TRUE)) +
      length(x) * (log(s1) + log(s2))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(b1 = c(0, 0), b2 = c(0, 0),
                              ls1 = 0, ls2 = 0, tr = 0),
                         silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)

  C <- rescor_matrix(fit)
  # the RE-free fit absorbs the simulated group effects into the
  # residuals, diluting the residual correlation to its marginal value
  marg <- (0.5 * 0.8 * 0.5) / sqrt((0.8^2 + 0.7^2) * (0.5^2 + 0.4^2))
  expect_lt(abs(C[1, 2] - marg), 0.12)
})

test_that("|ID| correlates random effects across responses", {
  dd <- sim_mv_data(seed = 63, re_cor = 0.7, n_g = 40, n = 800)
  fit <- suppressWarnings(
    frm(mvbf(bf(y1 ~ x + (1 | p | g)) + gaussian(),
             bf(y2 ~ x + (1 | p | g)) + gaussian()), data = dd)
  )
  expect_length(fit$frame$re_blocks, 1)
  bk <- fit$frame$re_blocks[[1]]
  expect_identical(bk$dim, 2L)
  V <- VarCorr(fit)[[1]]
  rho <- stats::cov2cor(V)[1, 2]
  expect_lt(abs(rho - 0.7), 0.25)

  # hand-rolled joint reference
  y1 <- dd$y1; y2 <- dd$y2; x <- dd$x; gi <- as.integer(dd$g)
  n_g <- nlevels(dd$g)
  nll_ref <- function(p) {
    "c" <- RTMB::ADoverload("c")
    sd_u <- exp(p$lsd)
    rho <- p$tr / sqrt(1 + p$tr^2)
    C <- RTMB::matrix(c(1, rho, rho, 1), 2, 2)
    S <- C * (RTMB::matrix(sd_u, ncol = 1) %*% RTMB::matrix(sd_u, nrow = 1))
    U <- RTMB::matrix(p$u, n_g, 2)
    nll <- -sum(RTMB::dmvnorm(U, 0, S, log = TRUE))
    mu1 <- p$b1[1] + p$b1[2] * x + U[gi, 1]
    mu2 <- p$b2[1] + p$b2[2] * x + U[gi, 2]
    nll - sum(RTMB::dnorm(y1, mu1, exp(p$ls1), log = TRUE)) -
      sum(RTMB::dnorm(y2, mu2, exp(p$ls2), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(b1 = c(0, 0), b2 = c(0, 0), ls1 = 0, ls2 = 0,
                              lsd = c(0, 0), tr = 0,
                              u = numeric(n_g * 2)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-5)
})

test_that("multivariate predict targets responses; post-processing guards", {
  dd <- sim_mv_data()
  mv <- frm(mvbf(bf(y1 ~ x + (1 | g)) + gaussian(),
                 bf(y2 ~ x + (1 | g)) + gaussian()), data = dd)
  p1 <- predict(mv, resp = "y1")
  p2 <- predict(mv, resp = "y2")
  expect_length(p1, nrow(dd))
  expect_false(isTRUE(all.equal(p1, p2)))
  expect_equal(predict(mv, newdata = dd, resp = "y2"), p2,
               tolerance = 1e-8)
  expect_error(predict(mv, resp = "zzz"), "Unknown response")
  expect_error(fitted(mv), "multivariate")
  expect_error(simulate(mv), "multivariate")
})

test_that("rescor validation", {
  dd <- sim_mv_data()
  dd$yc <- rpois(nrow(dd), 3)
  expect_error(
    frm(mvbf(bf(y1 ~ x) + gaussian(), bf(yc ~ x) + poisson(),
             rescor = TRUE), data = dd, dry_run = "spec"),
    "gaussian")
  dd$w <- runif(nrow(dd))
  expect_error(
    frm(mvbf(bf(y1 | weights(w) ~ x) + gaussian(),
             bf(y2 ~ x) + gaussian(), rescor = TRUE), data = dd),
    "weights")
})
