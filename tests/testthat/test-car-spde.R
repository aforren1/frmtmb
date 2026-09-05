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
  # (hard-constrained) likelihood type = "esicar" fits is four orders
  # below the parameter's own standard error; con_sd walks it down
  # quadratically
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

test_that("esicar matches a hand-rolled HARD-constrained ICAR ML", {
  # brms's own esicar parameterization as the reference: Nloc - 1 free
  # values and the last minus their sum, so cov(field) = sdcar^2 A Vz A'
  # with A the constraint basis. frmtmb reaches the same model from the
  # other side, keeping Nloc coefficients, centering them in
  # expand_b() and leaving the component mean inert. The two have to
  # agree to optimizer noise, because the Laplace approximation is
  # exact here.
  s <- car_lattice_data(42)
  W <- s$W
  A <- rbind(diag(s$n - 1), -1)
  Vz <- solve(t(A) %*% s$L %*% A)
  hard <- marginal_ml(s$d$y, s$X, s$Z,
                      function(p) exp(2 * p[1]) * (A %*% Vz %*% t(A)),
                      c(0, 0))
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
             data = s$d)
  expect_lt(abs(as.numeric(logLik(fit)) - hard$logLik), 1e-7)
  expect_vector_equal(c(fit$estimates$theta, fit$estimates$betad),
                      hard$par, tol = 1e-5)
  # EXACTLY zero, not merely small: icar's residual on this fit is 1e-9
  expect_lt(abs(sum(ranef(fit)[[1]])), 1e-12)
  expect_equal(unname(VarCorr(fit)[[1]][1, 1]),
               exp(2 * fit$estimates$theta[1]))
})

test_that("esicar and icar are different models, by the predicted amount", {
  # They selected one density through 0.51.0. icar's soft constraint is
  # the same model as an extra random intercept of sd con_sd * sdcar,
  # so its likelihood sits BELOW esicar's by the design note's bias,
  # 4.0e-4 at the 1e-3 default on this lattice, and walks onto it
  # quadratically as con_sd shrinks. esicar does not move at all,
  # because con_sd only scales a coordinate its predictor never sees.
  s <- car_lattice_data(42)
  W <- s$W
  ll <- function(ty, cs) {
    as.numeric(logLik(frm(bf(y ~ x + car(W, gr = loc, type = ty,
                                         con_sd = cs)) + gaussian(),
                          data = s$d)))
  }
  gap <- ll("esicar", 1e-3) - ll("icar", 1e-3)
  expect_gt(gap, 1e-4)
  expect_lt(gap, 1e-3)
  # the invariance IS the constraint being exact rather than tight
  es <- vapply(c(1e-2, 1e-3, 1e-4), function(cs) ll("esicar", cs), 0)
  ic <- vapply(c(1e-2, 1e-3, 1e-4), function(cs) ll("icar", cs), 0)
  expect_lt(diff(range(es)), 1e-8)
  expect_gt(diff(range(ic)), 1e-2)
  # and each decade of con_sd buys icar two more digits toward esicar
  expect_lt(abs(es[[3]] - ic[[3]]), abs(es[[2]] - ic[[2]]) / 50)
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

test_that("the esicar post-fit surface reads the FULL field", {
  # esicar keeps Nloc coefficients and centers them on the way to the
  # predictor, so everything downstream sees one value per location,
  # exactly as icar does. What is measured here is that they all read
  # the CENTERED vector: a path that read the raw parameters instead
  # would be off by the inert component mean.
  s <- car_lattice_data(3)
  W <- s$W
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
             data = s$d)
  re <- ranef(fit)[[1]]
  expect_equal(dim(re), c(s$n, 1L))
  expect_equal(rownames(re), rownames(s$W))
  expect_lt(abs(sum(re)), 1e-12)
  expect_equal(colnames(VarCorr(fit)[[1]]), "sd(car)")
  # in-sample and newdata prediction see the same design, standard
  # errors included
  rows <- c(1L, 10L, nrow(s$d))
  p_in <- predict(fit, se.fit = TRUE)
  p_nd <- predict(fit, newdata = s$d[rows, ], se.fit = TRUE)
  expect_vector_equal(p_nd$fit, p_in$fit[rows], tol = 1e-10)
  expect_vector_equal(p_nd$se.fit, p_in$se.fit[rows], tol = 1e-8)
  # a draw is Nloc long and its field is on the constraint again
  set.seed(1)
  bdraw <- frmtmb:::draw_b(fit)
  expect_length(bdraw, s$n)
  cdraw <- frmtmb:::expand_b(fit$frame, bdraw, fit$estimates$theta)
  expect_lt(abs(sum(cdraw)), 1e-12)
  set.seed(1)
  expect_equal(dim(simulate(fit, nsim = 2)), c(nrow(s$d), 2L))
  # the importance correction refuses every car type by name already,
  # and esicar keeps the level-major layout that refusal predates
  expect_error(
    frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
        data = s$d, importance = 32),
    "cannot correct the 'car' structure")
})

test_that("esicar centers its own block and leaves the others alone", {
  # expand_b() now has two exceptions, and a model carrying an esicar
  # block AND an ordinary one is where a mix-up would show: the CAR
  # field must come out on the constraint and the iid intercepts must
  # come out untouched. The reference is the marginal ML of the sum of
  # both variance components, with the CAR one at its EXACT constrained
  # covariance pinv(L).
  set.seed(11)
  W <- lattice_W(4, 4)
  n <- nrow(W)
  L <- diag(rowSums(W)) - W
  ev <- eigen(L, symmetric = TRUE)
  pos <- ev$values > 1e-8 * max(ev$values)
  Lp <- ev$vectors[, pos] %*% diag(1 / ev$values[pos]) %*%
    t(ev$vectors[, pos])
  phi <- drop(crossprod(chol(Lp + diag(1e-10, n)), stats::rnorm(n)))
  phi <- phi - mean(phi)
  loc <- factor(rep(rownames(W), each = 8), levels = rownames(W))
  g2 <- factor(rep(paste0("g", 1:8), length.out = length(loc)))
  u <- stats::rnorm(8, 0, 0.6)
  d <- data.frame(loc = loc, g2 = g2, x = stats::rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] + u[as.integer(d$g2)] +
    stats::rnorm(nrow(d), 0, 0.5)
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar") + (1 | g2)) +
               gaussian(), data = d)
  Zl <- stats::model.matrix(~ loc - 1, d)
  Zg <- stats::model.matrix(~ g2 - 1, d)
  ref <- marginal_ml(d$y, stats::model.matrix(~x, d),
                     cbind(Zl, Zg),
                     function(p) {
                       as.matrix(Matrix::bdiag(exp(2 * p[1]) * Lp,
                                               exp(2 * p[2]) * diag(8)))
                     },
                     c(0, log(0.6), log(0.5)))
  expect_lt(abs(as.numeric(logLik(fit)) - ref$logLik), 1e-7)
  re <- ranef(fit)
  car_re <- re[[grep("car", names(re))]]
  iid_re <- re[[grep("g2", names(re), fixed = TRUE)]]
  expect_equal(nrow(car_re), n)
  expect_equal(nrow(iid_re), 8L)
  # the CAR block is on the constraint and the iid block is not, which
  # is the whole point of centering ONE block
  expect_lt(abs(sum(car_re)), 1e-12)
  expect_gt(abs(sum(iid_re)), 1e-9)
})

test_that("esicar constrains a disconnected graph PER COMPONENT", {
  # brms constrains the global sum only, which leaves an intrinsic
  # field with two components improper. Constraining each component is
  # the same model whenever there is one component and a proper one
  # when there is not.
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
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
             data = d)
  # the reference is the Moore-Penrose inverse of the Laplacian, which
  # IS the covariance of the field constrained on both null directions
  ev <- eigen(L, symmetric = TRUE)
  pos <- ev$values > 1e-8 * max(ev$values)
  Lp <- ev$vectors[, pos] %*% diag(1 / ev$values[pos]) %*%
    t(ev$vectors[, pos])
  ref <- marginal_ml(d$y, stats::model.matrix(~x, d),
                     stats::model.matrix(~ loc - 1, d),
                     function(p) exp(2 * p[1]) * Lp, c(0, 0))
  expect_lt(abs(as.numeric(logLik(fit)) - ref$logLik), 1e-7)
  re <- ranef(fit)[[1]][, 1]
  expect_lt(max(abs(c(sum(re[seq_len(n1)]),
                      sum(re[n1 + seq_len(n1)])))), 1e-12)
})

test_that("the objective expands esicar whatever the frame's flag says", {
  # frame_needs_expand() derives the answer from the BLOCKS instead of
  # trusting frame$has_expand. `[[` does not partially match, so a
  # frame serialized before that field existed hands back NULL, and a
  # gate that read NULL as FALSE would evaluate the block WITHOUT
  # centering while its density still carried the constrained
  # normalizer: two different models and no error.
  #
  # It is invisible at the mode, because the inert coordinate is zero
  # there, so the check has to be made OFF the mode. That is the regime
  # imp_frozen_proposal() and cluster_scores_at() work in, and both
  # rebuild the objective from a STORED fit$frame.
  s <- car_lattice_data(42)
  W <- s$W
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
             data = s$d)
  pl <- fit$obj$env$parList(fit$obj$env$last.par.best)
  ref <- frmtmb:::build_objective(fit$frame)
  # a 0.51.0 frame: the field is absent, not FALSE
  fr_old <- fit$frame
  fr_old[["has_expand"]] <- NULL
  expect_null(fr_old[["has_expand"]])
  # and the belt-and-braces case, both flags actively FALSE
  fr_off <- fit$frame
  fr_off[["has_expand"]] <- FALSE
  fr_off[["has_rr"]] <- FALSE
  pl_off <- pl
  pl_off$b <- pl$b + 0.37
  # the shift has to be big enough that a missed centering would show
  expect_gt(ref(pl_off) - ref(pl), 1e3)
  for (frx in list(fr_old, fr_off)) {
    g <- frmtmb:::build_objective(frx)
    expect_equal(g(pl), ref(pl))
    expect_equal(g(pl_off), ref(pl_off))
  }
})

test_that("con_sd leaves the esicar fit alone but not its standard errors", {
  # con_sd cannot move the likelihood: the coordinate it scales enters
  # no linear predictor. It DOES reach the delta method, because
  # lp_delta_A() pairs the Z columns with b through dc/db = I rather
  # than the centering projection, so the coordinate's variance lands
  # in every standard error as exactly con_sd^2.
  #
  # Pinned here with the number, not the direction, so that the exact
  # Jacobian (R/predict.R) has something to flip when it lands.
  s <- car_lattice_data(42)
  W <- s$W
  fit <- function(cs) {
    frm(bf(y ~ x + car(W, gr = loc, type = "esicar", con_sd = cs)) +
          gaussian(), data = s$d)
  }
  f3 <- fit(1e-3)
  f4 <- fit(1e-4)
  # what con_sd does NOT touch
  expect_lt(abs(as.numeric(logLik(f3)) - as.numeric(logLik(f4))), 1e-8)
  expect_vector_equal(f3$estimates$theta, f4$estimates$theta, tol = 1e-8)
  expect_equal(unname(VarCorr(f3)[[1]][1, 1]),
               exp(2 * f3$estimates$theta[1]))
  # what it does: the variance excess is EXACTLY the difference of the
  # two con_sd squared, which identifies the leak rather than merely
  # bounding it
  se3 <- predict(f3, se.fit = TRUE)$se.fit
  se4 <- predict(f4, se.fit = TRUE)$se.fit
  expect_lt(max(abs((se3^2 - se4^2) - (1e-3^2 - 1e-4^2))), 1e-10)
  # 1.3e-5 relative in the standard error at the default
  rel <- max(se3 / sqrt(se3^2 - 1e-3^2) - 1)
  expect_gt(rel, 1e-5)
  expect_lt(rel, 2e-5)
  # ranef's conditional SDs leak the same way, and move with con_sd
  sd3 <- attr(ranef(f3, condVar = TRUE)[[1]], "condSD")
  sd4 <- attr(ranef(f4, condVar = TRUE)[[1]], "condSD")
  expect_gt(max(sd3 - sd4), 0)
  expect_lt(max(sd3 / sd4 - 1), 1e-4)
})

test_that("esicar handles a SINGLETON component", {
  # A component of one level contributes zero free dimensions: the
  # centering b_i - b_i is identically zero, and n - c counts it. The
  # reference is the Moore-Penrose inverse of the Laplacian, which is
  # the constrained covariance on all three components at once.
  n <- 9L
  W <- matrix(0, n, n)
  for (i in c(1L, 2L, 3L, 5L, 6L, 7L)) {
    W[i, i + 1L] <- 1
    W[i + 1L, i] <- 1
  }
  lv <- paste0("L", seq_len(n))
  dimnames(W) <- list(lv, lv)
  set.seed(5)
  L <- diag(rowSums(W)) - W
  ev <- eigen(L, symmetric = TRUE)
  pos <- ev$values > 1e-8 * max(ev$values)
  expect_equal(sum(pos), n - 3L)
  Lp <- ev$vectors[, pos] %*% diag(1 / ev$values[pos]) %*%
    t(ev$vectors[, pos])
  phi <- drop(crossprod(chol(Lp + diag(1e-10, n)), stats::rnorm(n)))
  loc <- factor(rep(lv, each = 8), levels = lv)
  d <- data.frame(loc = loc, x = stats::rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] +
    stats::rnorm(nrow(d), 0, 0.4)
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
             data = d)
  a <- fit$frame$re_blocks[[1]]$aux_car
  expect_equal(a$n_comp, 3L)
  expect_equal(as.numeric(a$nj), c(4, 4, 1))
  re <- ranef(fit)[[1]][, 1]
  # the singleton is not merely small, it is the number zero
  expect_identical(re[[9]], 0)
  expect_lt(max(abs(c(sum(re[1:4]), sum(re[5:8])))), 1e-12)
  ref <- marginal_ml(d$y, stats::model.matrix(~x, d),
                     stats::model.matrix(~ loc - 1, d),
                     function(p) exp(2 * p[1]) * Lp, c(0, log(0.4)))
  expect_lt(abs(as.numeric(logLik(fit)) - ref$logLik), 1e-7)
  # icar cannot do either: its singleton is free and its likelihood
  # sits below the constrained one
  fi <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")) + gaussian(),
            data = d)
  expect_gt(abs(ranef(fi)[[1]][, 1][[9]]), 1e-9)
  expect_lt(as.numeric(logLik(fi)), as.numeric(logLik(fit)))
  # escar still refuses a zero-degree location by name; the intrinsic
  # types do not need to
  expect_error(frm(bf(y ~ x + car(W, gr = loc, type = "escar")) +
                     gaussian(), data = d),
               "at least one")
})

test_that("esicar's con_sd invariance is not a gaussian accident", {
  # The factorization does not depend on the family: the component
  # means are decoupled from the data term whatever it is. Under
  # poisson the Laplace approximation is no longer exact, and the
  # invariance still holds, which says it comes from the model and not
  # from the gaussian's exactness.
  lattice <- lattice_W(4, 4)
  n <- nrow(lattice)
  set.seed(21)
  K <- diag(rowSums(lattice)) - lattice + matrix(1 / (1e-3 * n)^2, n, n)
  phi <- 0.6 * drop(crossprod(chol(solve(K)), stats::rnorm(n)))
  loc <- factor(rep(rownames(lattice), each = 6),
                levels = rownames(lattice))
  d <- data.frame(loc = loc, x = stats::rnorm(length(loc)))
  d$y <- stats::rpois(nrow(d),
                      exp(1 + 0.3 * d$x + phi[as.integer(d$loc)]))
  W <- lattice
  ll <- function(ty, cs) {
    as.numeric(logLik(frm(bf(y ~ x + car(W, gr = loc, type = ty,
                                         con_sd = cs)) + poisson(),
                          data = d)))
  }
  es <- c(ll("esicar", 1e-2), ll("esicar", 1e-4))
  ic <- c(ll("icar", 1e-2), ll("icar", 1e-4))
  expect_lt(diff(range(es)), 1e-8)
  expect_gt(diff(range(ic)), 1e-4)
  # and the constraint is still exact off the gaussian
  fp <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + poisson(),
            data = d)
  expect_lt(abs(sum(ranef(fp)[[1]])), 1e-12)
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
