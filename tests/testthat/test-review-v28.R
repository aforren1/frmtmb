# Heavy reference-validation file: full runs happen locally and in CI
# (NOT_CRAN=true). CRAN only needs to see the package work, not the
# validation program re-derived.
skip_on_cran()

# Confirmed defects from the v0.28 wave review.
#
# The one that changed answers is A1: spde() derived its mesh node
# ordering from the grouping variable's labels, which the finite-element
# matrices know nothing about, so integer node ids were read
# lexicographically and the field was fitted against a permuted
# precision. The rest are diagnosis and validation surfaces that either
# over-claimed, under-checked, or repeated themselves.

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

v28_marginal_ml <- function(y, X, Z, Sigma_fn, start) {
  nll <- function(p) {
    if (any(!is.finite(p)) || any(abs(p) > 25)) return(1e10)
    V <- tryCatch(Z %*% Sigma_fn(p) %*% t(Z) +
                    exp(2 * p[length(p)]) * diag(length(y)),
                  error = function(e) NULL)
    if (is.null(V)) return(1e10)
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

# Run expr with mc.cores set, collecting every warning it raises: an
# rstan run emits its own diagnostics, and this file must not leave
# stray warnings in the suite summary.
v28_warnings_under_cores <- function(cores, expr) {
  old <- options(mc.cores = cores)
  on.exit(options(old), add = TRUE)
  msgs <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = msgs)
}

# --------------------------------------------------- A1 spde node order

# The alpha = 2 finite-element triple of a regular 1-D linear mesh
# (lumped mass, stiffness, and M2 = G C^-1 G), the same construction the
# main spde test uses.
chain_fem <- function(nn, h = 0.5) {
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
  list(c0 = C0, g1 = G, g2 = G %*% solve(C0) %*% G)
}

# 12 nodes, so the lexicographic reading ("1", "10", "11", "12", "2",
# ...) differs from the numeric one and the permutation bites.
spde_chain_data <- function(seed = 5, nn = 12, per = 6) {
  fem <- chain_fem(nn)
  Qf <- function(lt, lk) {
    k2 <- exp(2 * lk)
    exp(2 * lt) * (k2 * k2 * fem$c0 + 2 * k2 * fem$g1 + fem$g2)
  }
  set.seed(seed)
  u <- drop(crossprod(chol(solve(Qf(0, log(0.8)))), stats::rnorm(nn)))
  node <- rep(seq_len(nn), each = per)
  d <- data.frame(node = node, x = stats::rnorm(nn * per))
  d$y <- 0.3 + 0.4 * d$x + u[node] + stats::rnorm(nrow(d), 0, 0.3)
  list(d = d, fem = fem, Qf = Qf, nn = nn)
}

test_that("every accepted spde node spelling is the same fit", {
  s <- spde_chain_data()
  fem <- s$fem
  d <- s$d
  # the direct marginal ML of the SAME model: the Laplace approximation
  # is exact for a gaussian response, so this is the answer, not a
  # tolerance
  ref <- v28_marginal_ml(d$y, stats::model.matrix(~x, d),
                     stats::model.matrix(~ factor(node) - 1, d),
                     function(p) solve(s$Qf(p[1], p[2])), c(0, 0, 0))
  spellings <- list(
    integer = d,
    double = transform(d, node = as.numeric(node)),
    factor = transform(d, node = factor(node)),
    character = transform(d, node = as.character(node)),
    # the spelling that used to permute the mesh: lexicographic levels
    factor_chr = transform(d, node = factor(as.character(node)))
  )
  for (nm in names(spellings)) {
    f <- frm(bf(y ~ x + spde(fem, gr = node)) + gaussian(),
             data = spellings[[nm]])
    expect_lt(abs(as.numeric(logLik(f)) - ref$logLik), 1e-7)
    expect_vector_equal(c(f$estimates$theta, f$estimates$betad),
                        ref$par, tol = 1e-4)
  }
})

test_that("spde refuses labels it cannot read as mesh rows", {
  s <- spde_chain_data(nn = 6, per = 5)
  fem <- s$fem
  d <- s$d
  # the mesh matrices carry no dimnames, so a label ordering is a guess
  # and the contract says so instead of guessing
  d_lab <- transform(d, node = paste0("L", node))
  expect_error(frm(bf(y ~ x + spde(fem, gr = node)) + gaussian(),
                   data = d_lab), "ROW NUMBER")
  d_f <- transform(d, node = factor(paste0("L", node)))
  expect_error(frm(bf(y ~ x + spde(fem, gr = node)) + gaussian(),
                   data = d_f), "whole-number node indices")
  # fractional indices are not rows either
  d_frac <- transform(d, node = node + 0.5)
  expect_error(frm(bf(y ~ x + spde(fem, gr = node)) + gaussian(),
                   data = d_frac), "whole-number node indices")
  # out of range, either end
  d_hi <- transform(d, node = node + 1L)
  expect_error(frm(bf(y ~ x + spde(fem, gr = node)) + gaussian(),
                   data = d_hi), "outside the mesh")
  d_lo <- transform(d, node = node - 1L)
  expect_error(frm(bf(y ~ x + spde(fem, gr = node)) + gaussian(),
                   data = d_lo), "outside the mesh")
})

test_that("an spde mesh node the data never visits keeps its column", {
  s <- spde_chain_data()
  d <- s$d[s$d$node != 7L, ]
  f <- frm(bf(y ~ x + spde(s$fem, gr = node)) + gaussian(), data = d)
  # the block is the MESH, not the observed levels: node 7 is still
  # there, at its prior
  re <- ranef(f)[[1]]
  expect_equal(nrow(re), s$nn)
  expect_equal(rownames(re), as.character(seq_len(s$nn)))
  expect_equal(ngrps(f)[[1]], s$nn)
  rows <- c(1L, 20L, nrow(d))
  expect_vector_equal(predict(f, newdata = d[rows, ]),
                      predict(f)[rows], tol = 1e-10)
})

test_that("the spde finite-element matrices must agree on the mesh size", {
  d <- data.frame(y = stats::rnorm(20), node = rep(1:5, 4))
  ok <- list(c0 = diag(5), g1 = diag(5), g2 = diag(5))
  bad <- ok
  bad$g1 <- diag(4)
  expect_error(frm(bf(y ~ spde(bad, gr = node)) + gaussian(), data = d),
               "mesh nodes")
  sq <- ok
  sq$c0 <- matrix(1, 5, 4)
  expect_error(frm(bf(y ~ spde(sq, gr = node)) + gaussian(), data = d),
               "square matrix")
})

# ------------------------------------------------- A2 call-valued gr

test_that("car() and spde() accept a call-valued gr", {
  s <- spde_chain_data(nn = 8, per = 6)
  f_call <- frm(bf(y ~ x + spde(s$fem, gr = factor(node))) + gaussian(),
                data = s$d)
  f_col <- frm(bf(y ~ x + spde(s$fem, gr = node)) + gaussian(),
               data = s$d)
  expect_equal(as.numeric(logLik(f_call)), as.numeric(logLik(f_col)))

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

# ------------------------------------------- A3 cores vs options(mc.cores)


# --------------------------------------------- A4 fit_error_context

test_that("only numerical failures get the undefined-likelihood cause", {
  set.seed(8)
  dd <- data.frame(x = stats::rnorm(40))
  dd$y <- 1 + dd$x + stats::rnorm(40)
  # an unknown optimizer is a typo, not an unbounded likelihood
  e1 <- tryCatch(frm(bf(y ~ x) + gaussian(), data = dd,
                     control = frmtmb_control(optimizer = "nope")),
                 error = function(e) conditionMessage(e))
  expect_match(e1, "Unknown optimizer 'nope'")
  expect_false(grepl("likelihood was undefined", e1))
  # neither is a custom optimizer that broke its return contract
  e2 <- tryCatch(frm(bf(y ~ x) + gaussian(), data = dd,
                     control = frmtmb_control(
                       optimizer = function(par, fn, gr, lower, upper,
                                            control) list(par = par))),
                 error = function(e) conditionMessage(e))
  expect_match(e2, "must return par, objective, and convergence")
  expect_false(grepl("likelihood was undefined", e2))
  # the model is still named, so the user knows which fit died
  expect_match(e1, "raised while fitting")

  # a genuinely undefined objective keeps the diagnosis and the remedies
  spec <- frm(bf(y ~ x) + gaussian(), data = dd, dry_run = "spec")
  e3 <- tryCatch(
    frmtmb:::fit_error_context(spec, NULL, FALSE, frmtmb_control(), NULL,
                               NULL, stop("NA/NaN function evaluation")),
    error = function(e) conditionMessage(e))
  expect_match(e3, "likelihood was undefined or unbounded")
  expect_match(e3, "NA/NaN function evaluation")
  expect_match(e3, "frmtmb_control\\(optimizer")
})

# ------------------------------------------------ A5 escar at rho -> 1

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
  eig <- fit$frame$re_blocks[[1]]$aux_car$eigW
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

# ------------------------------------- A6 extreme-theta reads log sds only

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

test_that("a collapsed log sd is still reported, by name", {
  set.seed(12)
  dd <- data.frame(g = factor(rep(1:8, 6)), x = stats::rnorm(48))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(48, 0, 0.7)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  # force the collapsed mode the heuristic is about
  fit$estimates$theta[1] <- -12
  dg <- diagnose(fit, quiet = TRUE)
  expect_length(dg$extreme_theta, 1L)
  # the message names the confint() row, not a bare theta index
  expect_match(names(dg$extreme_theta), "1 | g", fixed = TRUE)
  expect_match(names(dg$extreme_theta), "(Intercept)", fixed = TRUE)
  out <- utils::capture.output(diagnose(fit))
  expect_true(any(grepl("Extreme covariance parameters", out)))
  expect_true(any(grepl("1 | g", out, fixed = TRUE)))
})

# --------------------------------------------- A7 car adjacency validation

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

# --------------------------------------- A8 one non-finite-cov warning

test_that("a degenerate covariance warns once per fit", {
  cache <- new.env(parent = emptyenv())
  expect_warning(frmtmb:::warn_nonfinite_cov(cache), "not finite")
  expect_silent(frmtmb:::warn_nonfinite_cov(cache))
  # a fresh fit (or a refit, which replaces the cache) says it again
  expect_warning(frmtmb:::warn_nonfinite_cov(new.env(parent = emptyenv())),
                 "not finite")
  # no cache, no memory: the helper still warns when called bare
  expect_warning(frmtmb:::warn_nonfinite_cov(), "not finite")
})

test_that("repeated vcov() calls on a degenerate fit warn once", {
  set.seed(13)
  dd <- data.frame(x = stats::rnorm(30))
  dd$y <- 1 + dd$x + stats::rnorm(30)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  # a NaN covariance is what sdreport returns on a singular Hessian;
  # inject it so the degradation path runs without a pathological model
  sdr <- sdr_of(fit)
  sdr$cov.fixed[] <- NaN
  fit$cache$sdr <- sdr
  expect_warning(vcov(fit), "not finite")
  expect_silent(vcov(fit))
  expect_silent(summary(fit))
})
