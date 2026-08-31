# v0.11: ou/toep covariance structures, Kronecker gr(cov=) slopes,
# scalar-on-function regression, Gauss-Kronrod quadrature, pluggable
# optimizers.

test_that("built-in and custom optimizers agree", {
  data(sleepstudy, package = "lme4")
  bfm <- bf(Reaction ~ Days + (Days | Subject)) + gaussian()
  f1 <- frm(bfm, sleepstudy)
  f2 <- suppressWarnings(
    frm(bfm, sleepstudy, control = frmtmb_control(optimizer = "optim"))
  )   # L-BFGS-B stops a touch earlier; loglik still matches to 1e-4
  my_opt <- function(par, fn, gr, lower, upper, control) {
    r <- stats::optim(par, fn, gr, method = "BFGS",
                      control = list(maxit = 1000, reltol = 1e-12))
    list(par = r$par, objective = r$value, convergence = r$convergence)
  }
  f3 <- frm(bfm, sleepstudy,
            control = frmtmb_control(optimizer = my_opt))
  expect_lt(abs(as.numeric(logLik(f1)) - as.numeric(logLik(f2))), 1e-4)
  expect_lt(abs(as.numeric(logLik(f1)) - as.numeric(logLik(f3))), 1e-4)
  expect_error(frm(bfm, sleepstudy,
                   control = frmtmb_control(optimizer = "bogus")),
               "Unknown optimizer")
})

test_that("ou matches glmmTMB on irregular times", {
  skip_if_not_installed("glmmTMB")
  set.seed(701)
  n_g <- 50
  times <- c(0, 0.3, 0.4, 1.1, 1.5, 2.7)
  n_t <- length(times)
  rate <- 1.2; sd_u <- 0.9
  Sig <- sd_u^2 * exp(-rate * abs(outer(times, times, "-")))
  U <- matrix(rnorm(n_g * n_t), n_g) %*% chol(Sig)
  dd <- data.frame(
    y = 1 + as.vector(t(U)) + rnorm(n_g * n_t, 0, 0.4),
    g = factor(rep(seq_len(n_g), each = n_t)),
    tim = num_factor(rep(times, n_g))
  )
  fit <- frm(bf(y ~ 1 + ou(tim + 0 | g)) + gaussian(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ 1 + ou(tim + 0 | g), data = dd,
                          REML = FALSE)
  expect_loglik_equal(fit, ref, tol = 1e-5)
  # marginal sd recovered
  V <- VarCorr(fit)[[1]]
  expect_lt(abs(sqrt(V[1, 1]) - sd_u), 0.25)
  # decay with distance
  expect_gt(V[1, 2] / V[1, 1], V[1, 6] / V[1, 1])
})

test_that("ou validation needs numeric-coded coordinates", {
  dd <- sim_ar1_data(n_g = 5)
  dd$tim <- factor(letters[as.integer(dd$tim)])
  expect_error(frm(bf(y ~ 1 + ou(tim + 0 | g)) + gaussian(), data = dd),
               "num_factor")
})

test_that("toep matches a hand-rolled reference (glmmTMB when it converges)", {
  dd <- sim_ar1_data(seed = 702, n_g = 80, n_t = 4, rho = 0.5)
  fit <- suppressWarnings(
    frm(bf(y ~ 1 + toep(tim + 0 | g)) + gaussian(), data = dd)
  )

  yv <- dd$y; gi <- as.integer(dd$g); ti <- as.integer(dd$tim)
  d <- 4; ng <- nlevels(dd$g)
  nll_ref <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    sdv <- exp(p$lsd)
    cvec <- rep(p$phi[1], d)
    cvec[1] <- 1
    for (k in 1:3) cvec[k + 1] <- p$phi[k] / sqrt(1 + p$phi[k]^2)
    M <- abs(outer(1:d, 1:d, "-")) + 1L
    C <- RTMB::matrix(cvec[as.vector(M)], d, d)
    S <- C * (RTMB::matrix(sdv, ncol = 1) %*% RTMB::matrix(sdv, nrow = 1))
    U <- p$u
    dim(U) <- c(d, ng)
    nll <- -sum(RTMB::dmvnorm(t(U), 0, S, log = TRUE))
    mu <- p$b0 + p$u[(gi - 1) * d + ti]
    nll - sum(RTMB::dnorm(yv, mu, exp(p$ls), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(b0 = 0, ls = 0, lsd = numeric(d),
                              phi = numeric(3), u = numeric(d * ng)),
                         random = "u", silent = TRUE)
  opt <- suppressWarnings(
    nlminb(obj$par, obj$fn, obj$gr,
           control = list(iter.max = 2000, eval.max = 2000))
  )
  # the naive reference passes through non-PD proposals and stops a hair
  # short; equality to 1e-3 still pins the model down
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-3)

  # glmmTMB's toep is convergence-fragile on this data (false
  # convergence); compare only when it produces a finite likelihood
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    ref <- suppressWarnings(
      glmmTMB::glmmTMB(y ~ 1 + toep(tim + 0 | g), data = dd,
                       REML = FALSE)
    )
    ll_ref <- suppressWarnings(as.numeric(logLik(ref)))
    if (is.finite(ll_ref)) {
      expect_lt(abs(as.numeric(logLik(fit)) - ll_ref), 1e-4)
    }
  }
})

test_that("gr(cov=) with correlated slopes matches a Kronecker reference", {
  set.seed(703)
  ng <- 20
  R <- crossprod(matrix(rnorm(ng * ng), ng)) / ng
  A <- stats::cov2cor(R)
  dimnames(A) <- list(as.character(seq_len(ng)),
                      as.character(seq_len(ng)))
  Sig <- matrix(c(0.8^2, 0.3 * 0.8 * 0.4, 0.3 * 0.8 * 0.4, 0.4^2), 2)
  bvec <- drop(crossprod(chol(kronecker(A, Sig)), rnorm(2 * ng)))
  B <- matrix(bvec, nrow = 2)     # coefs x levels, level-major
  n <- 800
  g <- factor(rep(seq_len(ng), each = n / ng))
  x <- rnorm(n)
  dd <- data.frame(
    y = 1 + 0.5 * x + B[1, as.integer(g)] + B[2, as.integer(g)] * x +
      rnorm(n, 0, 0.5),
    x = x, g = g
  )

  fit <- frm(bf(y ~ x + (x | gr(g, cov = A))) + gaussian(), data = dd)
  expect_identical(fit$frame$re_blocks[[1]]$dim, 2L)

  yv <- dd$y; xv <- dd$x; gi <- as.integer(dd$g)
  Ad <- unname(A)
  nll_ref <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    sdv <- exp(p$ls2)
    C <- frmtmb:::us_chol_cor(p$cor, 2)
    S <- C * (RTMB::matrix(sdv, ncol = 1) %*% RTMB::matrix(sdv, nrow = 1))
    # build kron(A, S) with data index maps
    D <- 2 * ng
    r <- seq_len(D)
    l1 <- (r - 1) %/% 2 + 1; c1 <- (r - 1) %% 2 + 1
    ia <- as.vector(outer(l1, l1, function(a, b) (b - 1) * ng + a))
    is <- as.vector(outer(c1, c1, function(a, b) (b - 1) * 2 + a))
    K <- RTMB::matrix(as.vector(Ad)[ia] * as.vector(S)[is], D, D)
    nll <- -sum(RTMB::dmvnorm(p$u, 0, K, log = TRUE))
    mu <- p$b[1] + p$b[2] * xv + p$u[(gi - 1) * 2 + 1] +
      p$u[(gi - 1) * 2 + 2] * xv
    nll - sum(RTMB::dnorm(yv, mu, exp(p$ls), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(b = c(0, 0), ls = 0, ls2 = c(0, 0),
                              cor = 0, u = numeric(2 * ng)),
                         random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 1000, eval.max = 1000))
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-5)
  expect_identical(dim(VarCorr(fit)[[1]]), c(2L, 2L))
})

test_that("scalar-on-function regression matches mgcv exactly", {
  set.seed(602)
  n <- 250; K <- 30
  tt <- seq(0, 1, length.out = K)
  Xc <- matrix(rnorm(n * K), n, K)
  for (i in seq_len(n)) Xc[i, ] <- cumsum(rnorm(K, 0, 0.3))
  y <- 1 + drop(Xc %*% sin(2 * pi * tt)) / K + rnorm(n, 0, 0.3)
  dd <- data.frame(y = y)
  dd$Tmat <- matrix(tt, n, K, byrow = TRUE)
  dd$Lmat <- Xc / K   # curve values times quadrature weights

  fit <- frm(bf(y ~ s(Tmat, by = Lmat)) + gaussian(), data = dd)
  ref <- mgcv::gam(y ~ s(Tmat, by = Lmat), data = dd, method = "ML")
  expect_lt(abs(as.numeric(logLik(fit)) - (-as.numeric(ref$gcv.ubre))),
            1e-4)
  expect_lt(max(abs(fitted(fit) - fitted(ref))), 1e-4)
})

test_that("quadrature = TRUE matches glmer(nAGQ = 25)", {
  skip_if_not_installed("lme4")
  set.seed(601)
  ng <- 100; per <- 4
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per)
  u <- rnorm(ng, 0, 1.2)
  dd <- data.frame(y = rbinom(ng * per, 1, plogis(-0.5 + 0.7 * x + u[g])),
                   x = x, g = g)

  fit_gk <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd,
                quadrature = TRUE)
  ref <- lme4::glmer(y ~ x + (1 | g), dd, family = binomial, nAGQ = 25)
  expect_lt(abs(as.numeric(logLik(fit_gk)) - as.numeric(logLik(ref))),
            1e-4)
  expect_vector_equal(fixef(fit_gk)$mu, lme4::fixef(ref), tol = 1e-3)
  vc <- sqrt(VarCorr(fit_gk)[[1]][1, 1])
  sd_ref <- as.numeric(attr(lme4::VarCorr(ref)$g, "stddev"))[1]
  expect_lt(abs(vc - sd_ref), 1e-3)

  # Laplace and quadrature genuinely differ here
  fit_lap <- frm(bf(y ~ x + (1 | g)) + binomial(), data = dd)
  expect_gt(abs(as.numeric(logLik(fit_gk)) -
                  as.numeric(logLik(fit_lap))), 0.5)

  # second AGHQ implementation: GLMMadaptive
  if (requireNamespace("GLMMadaptive", quietly = TRUE)) {
    ga <- GLMMadaptive::mixed_model(y ~ x, random = ~ 1 | g, data = dd,
                                    family = binomial(), nAGQ = 25)
    expect_lt(abs(as.numeric(logLik(fit_gk)) - as.numeric(logLik(ga))),
              1e-3)
    expect_vector_equal(fixef(fit_gk)$mu,
                        unname(GLMMadaptive::fixef(ga)), tol = 1e-2)
  }

  # guardrails
  expect_error(frm(bf(y ~ x + (x | g)) + binomial(), data = dd,
                   quadrature = TRUE), "scalar random")
})