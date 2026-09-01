# Spatial GMRF grammar: car() (brms spelling) and spde().
#
# The primary reference everywhere is a hand-rolled direct ML fit of the
# same marginal gaussian model - V = Z Sigma_b Z' + sigma^2 I with the
# SAME block covariance, profiled over beta. The Laplace approximation is
# exact for a gaussian response, so agreement has to be to optimizer
# noise, not to a tolerance.

# Rook-adjacency of an r x c lattice, with brms-style dimnames.
lattice_W <- function(r, c) {
  g <- expand.grid(r = seq_len(r), c = seq_len(c))
  n <- nrow(g)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) W[i, j] <- 1
    }
  }
  dimnames(W) <- list(paste0("L", seq_len(n)), paste0("L", seq_len(n)))
  W
}

# -2 log profile likelihood / 2 of y ~ N(X beta, Z Sigma(p) Z' + s^2 I).
marginal_ml <- function(y, X, Z, Sigma_fn, start) {
  nll <- function(p) {
    if (any(!is.finite(p)) || any(abs(p) > 25)) return(1e10)
    V <- Z %*% Sigma_fn(p) %*% t(Z) +
      exp(2 * p[length(p)]) * diag(length(y))
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e10)
    Xs <- backsolve(R, X, transpose = TRUE)
    bh <- solve(crossprod(Xs), crossprod(Xs, backsolve(R, y,
                                                       transpose = TRUE)))
    r <- y - X %*% bh
    0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(R))) +
             sum(backsolve(R, r, transpose = TRUE)^2))
  }
  o <- stats::optim(start, nll, method = "BFGS",
                    control = list(reltol = 1e-14, maxit = 500))
  list(logLik = -o$value, par = o$par)
}

car_lattice_data <- function(seed, r = 4, c = 4, per = 6, sd_car = 1.2,
                             sigma = 0.5, con_sd = 1e-3) {
  set.seed(seed)
  W <- lattice_W(r, c)
  n <- nrow(W)
  L <- diag(rowSums(W)) - W
  K <- L + matrix(1 / (con_sd * n)^2, n, n)
  phi <- sd_car * drop(crossprod(chol(solve(K)), stats::rnorm(n)))
  loc <- factor(rep(rownames(W), each = per), levels = rownames(W))
  d <- data.frame(loc = loc, x = stats::rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] +
    stats::rnorm(nrow(d), 0, sigma)
  list(d = d, W = W, L = L, K = K, n = n,
       X = stats::model.matrix(~x, d),
       Z = stats::model.matrix(~ loc - 1, d))
}

test_that("icar matches a hand-rolled constrained-ICAR ML", {
  s <- car_lattice_data(42)
  W <- s$W
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")) + gaussian(),
             data = s$d)
  Kinv <- solve(s$K)
  ref <- marginal_ml(s$d$y, s$X, s$Z,
                     function(p) exp(2 * p[1]) * Kinv, c(0, 0))
  expect_lt(abs(as.numeric(logLik(fit)) - ref$logLik), 1e-7)
  expect_vector_equal(c(fit$estimates$theta, fit$estimates$betad),
                      ref$par, tol = 1e-5)
  # the constraint does its job: the field sums to (near) zero. The
  # soft constraint pins the sum at an sd of con_sd * n * sdcar, and
  # the conditional mode lands far inside that
  expect_lt(abs(sum(ranef(fit)[[1]])), 1e-3)
  # and the fitted sdcar is the field scale, not a variance
  expect_equal(unname(VarCorr(fit)[[1]][1, 1]),
               exp(2 * fit$estimates$theta[1]))
})

test_that("the soft sum-to-zero constraint converges to the hard one", {
  # the default is brms's 1e-3, whose distance from the exact
  # (hard-constrained, brms esicar) likelihood is four orders below the
  # parameter's own standard error; con_sd walks it down quadratically
  s <- car_lattice_data(42)
  W <- s$W
  A <- rbind(diag(s$n - 1), -1)
  Vz <- solve(t(A) %*% s$L %*% A)
  hard <- marginal_ml(s$d$y, s$X, s$Z,
                      function(p) exp(2 * p[1]) * (A %*% Vz %*% t(A)),
                      c(0, 0))
  gap <- vapply(c(1e-3, 1e-4, 1e-5), function(cs) {
    f <- frm(bf(y ~ x + car(W, gr = loc, type = "icar", con_sd = cs)) +
               gaussian(), data = s$d)
    c(abs(as.numeric(logLik(f)) - hard$logLik),
      abs(exp(f$estimates$theta[1]) - exp(hard$par[1])) / exp(hard$par[1]))
  }, numeric(2))
  # the default is small, and each decade buys two more
  expect_lt(gap[1, 1], 1e-3)
  expect_lt(gap[2, 1], 1e-4)
  expect_lt(gap[1, 3], gap[1, 1] / 50)
  expect_lt(gap[2, 3], gap[2, 1] / 50)
  # the default IS the default
  f_def <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")) + gaussian(),
               data = s$d)
  f_1e3 <- frm(bf(y ~ x + car(W, gr = loc, type = "icar",
                              con_sd = 1e-3)) + gaussian(), data = s$d)
  expect_equal(as.numeric(logLik(f_def)), as.numeric(logLik(f_1e3)))
})

test_that("esicar selects the same density as icar", {
  s <- car_lattice_data(42)
  W <- s$W
  f1 <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")) + gaussian(),
            data = s$d)
  f2 <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
            data = s$d)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)))
  expect_vector_equal(f1$estimates$theta, f2$estimates$theta, tol = 1e-12)
})

test_that("escar matches a hand-rolled proper-CAR ML", {
  # data from a PROPER CAR, so the dependence parameter has an interior
  # optimum and both fits can be compared parameter by parameter
  set.seed(7)
  W <- lattice_W(6, 6)
  n <- nrow(W)
  Dg <- diag(rowSums(W))
  phi <- drop(crossprod(chol(solve(Dg - 0.7 * W)), stats::rnorm(n)))
  loc <- factor(rep(rownames(W), each = 8), levels = rownames(W))
  d <- data.frame(loc = loc, x = stats::rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] +
    stats::rnorm(nrow(d), 0, 0.5)
  s <- list(d = d, X = stats::model.matrix(~x, d),
            Z = stats::model.matrix(~ loc - 1, d))
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "escar")) + gaussian(),
             data = s$d)
  ref <- marginal_ml(s$d$y, s$X, s$Z, function(p) {
    exp(2 * p[1]) * solve(Dg - stats::plogis(p[2]) * W)
  }, c(0, 0, 0))
  expect_lt(abs(as.numeric(logLik(fit)) - ref$logLik), 1e-7)
  expect_vector_equal(c(fit$estimates$theta, fit$estimates$betad),
                      ref$par, tol = 1e-4)
  # the dependence parameter is reported on brms's (0, 1) scale
  cv <- confint_varcorr(fit)
  expect_true(all(c("sd(car)", "car") %in% cv$term))
  rho_row <- cv[cv$term == "car", ]
  expect_equal(rho_row$estimate, stats::plogis(fit$estimates$theta[2]))
  expect_true(rho_row$lwr > 0 && rho_row$upr < 1)
})

test_that("bym2 matches the hand-rolled scaled mixture", {
  set.seed(11)
  W <- lattice_W(6, 6)
  n <- nrow(W)
  L <- diag(rowSums(W)) - W
  Ki <- solve(L + matrix(1 / (1e-3 * n)^2, n, n))
  scl <- as.numeric(car_scale_factor(Matrix::Matrix(W, sparse = TRUE)))
  Sig <- 1 * ((1 - 0.6) * diag(n) + 0.6 / scl * Ki)
  phi <- drop(crossprod(chol(Sig), stats::rnorm(n)))
  loc <- factor(rep(rownames(W), each = 8), levels = rownames(W))
  d <- data.frame(loc = loc, x = stats::rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] +
    stats::rnorm(nrow(d), 0, 0.5)
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "bym2")) + gaussian(),
             data = d)
  ref <- marginal_ml(d$y, stats::model.matrix(~x, d),
                     stats::model.matrix(~ loc - 1, d), function(p) {
                       rho <- stats::plogis(p[2])
                       exp(2 * p[1]) * ((1 - rho) * diag(n) +
                                          rho / scl * Ki)
                     }, c(0, 0, 0))
  expect_lt(abs(as.numeric(logLik(fit)) - ref$logLik), 1e-7)
  expect_vector_equal(c(fit$estimates$theta, fit$estimates$betad),
                      ref$par, tol = 1e-4)
  expect_true("rhocar" %in% confint_varcorr(fit)$term)
})

test_that("the bym2 scaling factor is brms's", {
  # brms:::.car_scale, written against the edge list; ours reads the
  # adjacency matrix, so the two agree only if the construction does
  W <- lattice_W(3, 4)
  n <- nrow(W)
  ours <- as.numeric(car_scale_factor(Matrix::Matrix(W, sparse = TRUE)))
  Q <- Matrix::Diagonal(n, Matrix::rowSums(W)) - Matrix::Matrix(W)
  Qp <- Q + Matrix::Diagonal(n) * max(Matrix::diag(Q)) *
    sqrt(.Machine$double.eps)
  S <- Matrix::solve(Qp)
  A <- matrix(1, 1, n)
  Wc <- S %*% t(A)
  S <- S - Wc %*% solve(A %*% Wc) %*% Matrix::t(Wc)
  expect_equal(ours, exp(mean(log(Matrix::diag(S)))))
})

test_that("a disconnected graph gets the right rank correction", {
  # two 2 x 3 lattices with no edge between them: the intrinsic field
  # has a two-dimensional null space, so the constraint (and the
  # log-determinant) must be applied per component
  set.seed(99)
  W1 <- lattice_W(2, 3)
  n1 <- nrow(W1)
  W <- matrix(0, 2 * n1, 2 * n1)
  W[seq_len(n1), seq_len(n1)] <- W1
  W[n1 + seq_len(n1), n1 + seq_len(n1)] <- W1
  lv <- paste0("L", seq_len(2 * n1))
  dimnames(W) <- list(lv, lv)
  n <- 2 * n1
  L <- diag(rowSums(W)) - W
  S <- rbind(c(rep(1, n1), rep(0, n1)), c(rep(0, n1), rep(1, n1)))
  K <- L + t(S) %*% diag(rep(1 / (1e-3 * n1)^2, 2)) %*% S
  phi <- drop(crossprod(chol(solve(K)), stats::rnorm(n)))
  loc <- factor(rep(lv, each = 8), levels = lv)
  d <- data.frame(loc = loc, x = stats::rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] +
    stats::rnorm(nrow(d), 0, 0.4)
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")) + gaussian(),
             data = d)
  ref <- marginal_ml(d$y, stats::model.matrix(~x, d),
                     stats::model.matrix(~ loc - 1, d),
                     function(p) exp(2 * p[1]) * solve(K), c(0, 0))
  expect_lt(abs(as.numeric(logLik(fit)) - ref$logLik), 1e-7)
  # each component sums to zero on its own
  re <- ranef(fit)[[1]][, 1]
  expect_lt(max(abs(c(sum(re[seq_len(n1)]), sum(re[n1 + seq_len(n1)])))),
            1e-3)
})

test_that("car recovers its parameters on repeated lattices", {
  # single-realization spatial fits are noisy, so recovery is judged on
  # the mean over replicates
  set.seed(2024)
  W <- lattice_W(5, 5)
  n <- nrow(W)
  L <- diag(rowSums(W)) - W
  Ki <- solve(L + matrix(1 / (1e-3 * n)^2, n, n))
  Rk <- chol(Ki)
  sd_true <- 1.0
  est <- vapply(seq_len(15), function(i) {
    phi <- sd_true * drop(crossprod(Rk, stats::rnorm(n)))
    loc <- factor(rep(rownames(W), each = 8), levels = rownames(W))
    d <- data.frame(loc = loc)
    d$y <- 1 + phi[as.integer(d$loc)] + stats::rnorm(nrow(d), 0, 0.5)
    f <- frm(bf(y ~ 1 + car(W, gr = loc, type = "icar")) + gaussian(),
             data = d)
    exp(f$estimates$theta[1])
  }, numeric(1))
  expect_lt(abs(mean(est) - sd_true), 0.15)
})

test_that("bym2 recovers its mixing parameter on repeated lattices", {
  set.seed(4242)
  W <- lattice_W(5, 5)
  n <- nrow(W)
  L <- diag(rowSums(W)) - W
  Ki <- solve(L + matrix(1 / (1e-3 * n)^2, n, n))
  scl <- as.numeric(car_scale_factor(Matrix::Matrix(W, sparse = TRUE)))
  rho_true <- 0.7
  sd_true <- 1.0
  Rk <- chol(sd_true^2 * ((1 - rho_true) * diag(n) + rho_true / scl * Ki))
  est <- vapply(seq_len(15), function(i) {
    phi <- drop(crossprod(Rk, stats::rnorm(n)))
    loc <- factor(rep(rownames(W), each = 10), levels = rownames(W))
    d <- data.frame(loc = loc)
    d$y <- 1 + phi[as.integer(d$loc)] + stats::rnorm(nrow(d), 0, 0.4)
    f <- frm(bf(y ~ 1 + car(W, gr = loc, type = "bym2")) + gaussian(),
             data = d)
    c(exp(f$estimates$theta[1]), stats::plogis(f$estimates$theta[2]))
  }, numeric(2))
  expect_lt(abs(mean(est[1, ]) - sd_true), 0.15)
  expect_lt(abs(mean(est[2, ]) - rho_true), 0.2)
})

test_that("the car post-fit surface answers", {
  s <- car_lattice_data(3)
  W <- s$W
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "bym2")) + gaussian(),
             data = s$d)
  re <- ranef(fit)[[1]]
  expect_equal(dim(re), c(s$n, 1L))
  expect_equal(rownames(re), rownames(s$W))
  expect_equal(colnames(VarCorr(fit)[[1]]), "sd(car)")
  cv <- confint_varcorr(fit)
  expect_true(all(cv$lwr < cv$estimate & cv$estimate < cv$upr))
  # in-sample and newdata prediction see the same design
  rows <- c(1L, 10L, nrow(s$d))
  p_in <- predict(fit, se.fit = TRUE)
  p_nd <- predict(fit, newdata = s$d[rows, ], se.fit = TRUE)
  expect_vector_equal(p_nd$fit, p_in$fit[rows], tol = 1e-10)
  expect_vector_equal(p_nd$se.fit, p_in$se.fit[rows], tol = 1e-8)
  expect_equal(ngrps(fit)[["loc"]], s$n)
  set.seed(1)
  expect_equal(dim(simulate(fit, nsim = 2)), c(nrow(s$d), 2L))
  set.seed(1)
  expect_length(frmtmb:::draw_b(fit), s$n)
  # a location the fit never saw has no structure to borrow
  nd <- s$d[1:3, ]
  nd$loc <- factor("ZZ", levels = c(levels(s$d$loc), "ZZ"))
  expect_error(predict(fit, newdata = nd), "New levels")
  expect_silent(predict(fit, newdata = nd, allow_new_levels = TRUE))
})

test_that("car validates its adjacency matrix and its grammar", {
  s <- car_lattice_data(3, r = 2, c = 3, per = 4)
  d <- s$d
  W <- s$W
  bad <- W
  bad[1, 2] <- 5
  expect_error(frm(bf(y ~ car(bad, gr = loc)) + gaussian(), data = d),
               "symmetric")
  asym <- W
  dimnames(asym) <- NULL
  expect_error(frm(bf(y ~ car(asym, gr = loc)) + gaussian(), data = d),
               "dimnames")
  short <- W[-1, -1]
  expect_error(frm(bf(y ~ car(short, gr = loc)) + gaussian(), data = d),
               "no row for location")
  expect_error(frm(bf(y ~ car(W, gr = loc, type = "bym")) + gaussian(),
                   data = d), "type must be one of")
  expect_error(frm(bf(y ~ car(W, gr = loc, foo = 1)) + gaussian(),
                   data = d), "unknown argument")
  expect_error(frm(bf(y ~ (1 | car(W, gr = loc))) + gaussian(), data = d),
               "not a bar term")
  # brms's deprecated gr = NA default is refused by name
  expect_error(frm(bf(y ~ car(W, gr = NA)) + gaussian(), data = d),
               "gr must name a grouping variable")
  # weighted adjacency is binarized, brms-style, with a message
  wt <- W * 2
  expect_message(frm(bf(y ~ car(wt, gr = loc, type = "icar")) +
                       gaussian(), data = d), "non-zero values")
  # an isolated location has no proper CAR conditional
  iso <- W
  iso[1, ] <- 0
  iso[, 1] <- 0
  expect_error(frm(bf(y ~ car(iso, gr = loc, type = "escar")) +
                     gaussian(), data = d), "at least one neighbor")
})

test_that("spde matches a dense direct ML on a 1-D chain", {
  # the alpha = 2 finite-element triple of a regular 1-D linear mesh:
  # lumped mass C, stiffness G, and M2 = G C^-1 G. Same identities the
  # planar fmesher/INLA matrices satisfy, so the assembly is exercised
  # exactly as it would be on a real mesh.
  set.seed(5)
  nn <- 20
  h <- 0.5
  C0 <- diag(rep(h, nn))
  C0[1, 1] <- h / 2
  C0[nn, nn] <- h / 2
  G <- matrix(0, nn, nn)
  for (i in seq_len(nn - 1)) {
    G[i, i] <- G[i, i] + 1 / h
    G[i + 1, i + 1] <- G[i + 1, i + 1] + 1 / h
    G[i, i + 1] <- -1 / h
    G[i + 1, i] <- -1 / h
  }
  fem <- list(c0 = C0, g1 = G, g2 = G %*% solve(C0) %*% G)
  Qf <- function(lt, lk) {
    k2 <- exp(2 * lk)
    exp(2 * lt) * (k2 * k2 * C0 + 2 * k2 * G + fem$g2)
  }
  u <- drop(crossprod(chol(solve(Qf(0, log(0.8)))), stats::rnorm(nn)))
  node <- factor(rep(seq_len(nn), each = 6))
  d <- data.frame(node = node, x = stats::rnorm(nn * 6))
  d$y <- 0.3 + 0.4 * d$x + u[as.integer(d$node)] +
    stats::rnorm(nrow(d), 0, 0.3)
  fit <- frm(bf(y ~ x + spde(fem, gr = node)) + gaussian(), data = d)
  ref <- marginal_ml(d$y, stats::model.matrix(~x, d),
                     stats::model.matrix(~ node - 1, d),
                     function(p) solve(Qf(p[1], p[2])), c(0, 0, 0))
  expect_lt(abs(as.numeric(logLik(fit)) - ref$logLik), 1e-7)
  expect_vector_equal(c(fit$estimates$theta, fit$estimates$betad),
                      ref$par, tol = 1e-4)
  # the INLA spelling of the same matrices is the same model
  fem2 <- list(M0 = fem$c0, M1 = fem$g1, M2 = fem$g2)
  fit2 <- frm(bf(y ~ x + spde(fem2, gr = node)) + gaussian(), data = d)
  expect_equal(as.numeric(logLik(fit2)), as.numeric(logLik(fit)))
  # post-fit surface
  expect_equal(dim(ranef(fit)[[1]]), c(nn, 1L))
  expect_equal(sort(confint_varcorr(fit)$term),
               c("range(spde)", "sd(spde)"))
  expect_equal(colnames(VarCorr(fit)[[1]]), "sd(spde)")
  rows <- c(2L, 40L, 100L)
  expect_vector_equal(predict(fit, newdata = d[rows, ], se.fit = TRUE)$fit,
                      predict(fit, se.fit = TRUE)$fit[rows], tol = 1e-10)
  set.seed(1)
  expect_length(frmtmb:::draw_b(fit), nn)
})

test_that("spde validates its finite-element matrices", {
  d <- data.frame(y = stats::rnorm(20), node = factor(rep(1:5, 4)))
  ok <- list(c0 = diag(5), g1 = diag(5), g2 = diag(5))
  expect_error(frm(bf(y ~ spde(list(a = diag(5)), gr = node)) +
                     gaussian(), data = d), "M0, M1, M2")
  bad <- ok
  bad$g1 <- diag(4)
  expect_error(frm(bf(y ~ spde(bad, gr = node)) + gaussian(), data = d),
               "mesh nodes")
  expect_error(frm(bf(y ~ (1 | spde(ok, gr = node))) + gaussian(),
                   data = d), "not a bar term")
})

test_that("the spatial GMRF blocks are declared in the registry", {
  for (cs in c("car", "spde")) {
    expect_equal(frm_compat(cs, "gaussian")$status, "conditional",
                 info = cs)
    expect_equal(frm_compat(cs, "simulate")$status, "works", info = cs)
  }
  expect_match(frm_compat("car", "poisson")$note, "adjacency")
  expect_match(frm_compat("spde", "poisson")$note, "finite-element")
  expect_match(frm_compat("gr_prec", "gaussian")$note, "Kronecker")
})
