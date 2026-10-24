# car() blocks moved from frmtmb's tests/testthat/test-review-v28.R (the
# v28 review's car items A5 and A7 and its boundary check), unchanged.

# Test files see the package namespace but not each other, so the two
# references this file shares with test-car-spde.R are repeated here:
# the rook adjacency of a lattice, and the direct marginal ML of
# y ~ N(X beta, Z Sigma(p) Z' + s^2 I) profiled over beta. The Laplace
# approximation is exact for a gaussian response, so that ML value is
# the answer a fit has to reproduce, not a tolerance to sit inside.
v28_lattice_W <- function(r, c) {
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

test_that("escar stays finite as its dependence parameter reaches 1", {
  set.seed(7)
  W <- v28_lattice_W(4, 4)
  n <- nrow(W)
  loc <- factor(rep(rownames(W), each = 6), levels = rownames(W))
  dd <- data.frame(loc = loc, x = stats::rnorm(length(loc)))
  dd$y <- 1 + 0.5 * dd$x +
    stats::rnorm(n, 0, 0.8)[as.integer(dd$loc)] +
    stats::rnorm(nrow(dd), 0, 0.5)
  fit <- frm(bf(y ~ x + car(W, gr = loc, type = "escar")) + gaussian(),
             data = dd)
  # the normalized adjacency's spectrum is theoretically in [-1, 1] and
  # LAPACK returns it a few ulp outside; log(1 - rho * eig) with
  # rho = plogis(36) = 1 to double precision is then NaN
  eig <- fit$frame$re_blocks[[1]]$aux$eigW
  expect_true(all(eig >= -1 & eig <= 1))
  p <- fit$opt$par
  p[which(names(p) == "theta")[2]] <- 36
  expect_true(is.finite(fit$obj$fn(p)))
  expect_true(all(is.finite(fit$obj$gr(p))))
  # and the clamp costs nothing away from the boundary: the CAR
  # log-determinant it feeds is the exact one to floating point
  deg <- rowSums(W)
  isq <- diag(1 / sqrt(deg), nrow = n)
  raw <- eigen(isq %*% W %*% isq, symmetric = TRUE,
               only.values = TRUE)$values
  for (rho in c(0.1, 0.5, 0.9, 0.999)) {
    expect_equal(sum(log(1 - rho * eig)), sum(log(1 - rho * raw)),
                 tolerance = 1e-12)
  }
  # the exact determinant, straight from the precision matrix
  rho <- 0.9
  ldet <- as.numeric(Matrix::determinant(
    Matrix::Matrix(diag(deg) - rho * W), logarithm = TRUE)$modulus)
  expect_equal(sum(log(deg)) + sum(log(1 - rho * eig)), ldet,
               tolerance = 1e-10)
})

test_that("boundary parameters that are not log sds are not singular", {
  # a bym2 field with no independent part drives the mixing parameter to
  # its (0, 1) boundary, which is theta = +20 or more on the internal
  # scale - an ordinary optimum, not a collapsed variance
  set.seed(11)
  W <- v28_lattice_W(4, 4)
  n <- nrow(W)
  L <- diag(rowSums(W)) - W
  Ki <- solve(L + matrix(1 / (1e-3 * n)^2, n, n))
  scl <- as.numeric(car_scale_factor(Matrix::Matrix(W, sparse = TRUE)))
  phi <- drop(crossprod(chol(1.2^2 * Ki / scl), stats::rnorm(n)))
  loc <- factor(rep(rownames(W), each = 8), levels = rownames(W))
  db <- data.frame(loc = loc)
  db$y <- 1 + phi[as.integer(db$loc)] + stats::rnorm(nrow(db), 0, 0.3)
  fb <- frm(bf(y ~ 1 + car(W, gr = loc, type = "bym2")) + gaussian(),
            data = db)
  expect_gt(abs(fb$estimates$theta[2]), 8)
  expect_lt(abs(fb$estimates$theta[1]), 8)
  expect_length(diagnose(fb, quiet = TRUE)$extreme_theta, 0L)
  # the printed report is silent about it too
  out <- utils::capture.output(diagnose(fb))
  expect_false(any(grepl("Extreme covariance", out)))
  # frm_sample's mode-init warning reads the same components, so the
  # mixing parameter cannot trip it either
  idx <- log_sd_theta_index(fb)
  expect_equal(unname(idx), 1L)
  expect_match(names(idx), "sd\\(car\\)")
  expect_length(idx[abs(fb$estimates$theta[idx]) > 8], 0L)
})

test_that("car() validates its adjacency matrix's contents", {
  set.seed(3)
  W <- v28_lattice_W(2, 3)
  loc <- factor(rep(rownames(W), each = 5), levels = rownames(W))
  dd <- data.frame(loc = loc)
  dd$y <- stats::rnorm(nrow(dd))
  # a missing entry used to reach the user as "missing value where
  # TRUE/FALSE needed" from isSymmetric()
  Wna <- W
  Wna[2, 3] <- NA
  expect_error(frm(bf(y ~ car(Wna, gr = loc)) + gaussian(), data = dd),
               "missing entry")
  # a negative weight is a different model; brms's validate_car_matrix
  # refuses it rather than binarizing it into a neighbor
  Wneg <- W
  Wneg[1, 2] <- -1
  Wneg[2, 1] <- -1
  expect_error(frm(bf(y ~ car(Wneg, gr = loc)) + gaussian(), data = dd),
               "negative entry")
  # duplicate location names cannot be matched
  Wdup <- W
  dimnames(Wdup) <- list(c("a", "a", "b", "c", "d", "e"),
                         c("a", "a", "b", "c", "d", "e"))
  dd2 <- data.frame(loc = factor(rep(c("a", "b", "c", "d", "e"), each = 5)))
  dd2$y <- stats::rnorm(nrow(dd2))
  expect_error(frm(bf(y ~ car(Wdup, gr = loc)) + gaussian(), data = dd2),
               "more than once")
  # rownames and colnames that disagree name two different places
  Wm <- W
  colnames(Wm) <- rev(rownames(W))
  expect_error(frm(bf(y ~ car(Wm, gr = loc)) + gaussian(), data = dd),
               "same locations")
})

test_that("car() accepts a one-sided set of dimnames", {
  set.seed(4)
  W <- v28_lattice_W(2, 3)
  loc <- factor(rep(rownames(W), each = 6), levels = rownames(W))
  dd <- data.frame(loc = loc)
  dd$y <- 1 + stats::rnorm(6, 0, 0.6)[as.integer(dd$loc)] +
    stats::rnorm(nrow(dd), 0, 0.5)
  ll <- function(M) {
    as.numeric(logLik(frm(bf(y ~ car(M, gr = loc, type = "icar")) +
                            gaussian(), data = dd)))
  }
  Wr <- W
  colnames(Wr) <- NULL
  Wc <- W
  rownames(Wc) <- NULL
  expect_equal(ll(Wr), ll(W))
  expect_equal(ll(Wc), ll(W))
  Wnone <- W
  dimnames(Wnone) <- NULL
  expect_error(frm(bf(y ~ car(Wnone, gr = loc)) + gaussian(), data = dd),
               "dimnames")
})

test_that("car() accepts a call-valued gr", {
  set.seed(2)
  W <- matrix(0, 6, 6)
  for (i in 1:5) {
    W[i, i + 1] <- 1
    W[i + 1, i] <- 1
  }
  dimnames(W) <- list(paste0("L", 1:6), paste0("L", 1:6))
  dc <- data.frame(lab = rep(paste0("L", 1:6), each = 10))
  dc$y <- stats::rnorm(60, stats::rnorm(6, 0, 0.8)[
    as.integer(factor(dc$lab))], 0.5)
  dc$labf <- factor(dc$lab)
  c_call <- frm(bf(y ~ car(W, gr = factor(lab), type = "icar")) +
                  gaussian(), data = dc)
  c_col <- frm(bf(y ~ car(W, gr = labf, type = "icar")) + gaussian(),
               data = dc)
  expect_equal(as.numeric(logLik(c_call)), as.numeric(logLik(c_col)))
  expect_vector_equal(c_call$estimates$theta, c_col$estimates$theta,
                      tol = 1e-10)
})

