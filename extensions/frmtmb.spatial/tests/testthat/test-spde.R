# spde(), moved from frmtmb's tests/testthat/test-car-spde.R, test-review-v28.R,
# test-data2.R and test-tmb-examples.R with the test_that() blocks unchanged.
#
# The primary reference is a hand-rolled direct ML fit of the same marginal
# gaussian model, V = Z Sigma_b Z' + sigma^2 I with the SAME block
# covariance, profiled over beta. The Laplace approximation is exact for
# a gaussian response, so agreement has to be to optimizer noise.

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

# alpha = 2 finite-element triple of a regular 1-D linear mesh.
d2_chain_fem <- function(nn, h = 0.5) {
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
  expect_equal(colnames(varcorr_matrices(fit)[[1]]), "sd(spde)")
  rows <- c(2L, 40L, 100L)
  expect_vector_equal(frm_linpred(fit, newdata = d[rows, ], se.fit = TRUE)$fit,
                      frm_linpred(fit, se.fit = TRUE)$fit[rows], tol = 1e-10)
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
  # car() is frmtmb.spatial's, whose test-car.R checks its rows
  for (cs in c("spde")) {
    expect_equal(frm_compat(cs, "gaussian")$status, "conditional",
                 info = cs)
    expect_equal(frm_compat(cs, "simulate")$status, "works", info = cs)
  }
  expect_match(frm_compat("spde", "poisson")$note, "finite-element")
  expect_match(frm_compat("gr_prec", "gaussian")$note, "Kronecker")
})

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
  # an spde block is not a brms grouping factor, so ngrps() has no
  # entry for it since item 2.6f
  expect_null(ngrps(f))
  expect_equal(f$frame[["re_blocks"]][[1L]][["n_levels"]], s$nn)
  rows <- c(1L, 20L, nrow(d))
  expect_vector_equal(frm_linpred(f, newdata = d[rows, ]),
                      frm_linpred(f)[rows], tol = 1e-10)
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

test_that("spde() accepts a call-valued gr", {
  s <- spde_chain_data(nn = 8, per = 6)
  f_call <- frm(bf(y ~ x + spde(s$fem, gr = factor(node))) + gaussian(),
                data = s$d)
  f_col <- frm(bf(y ~ x + spde(s$fem, gr = node)) + gaussian(),
               data = s$d)
  expect_equal(as.numeric(logLik(f_call)), as.numeric(logLik(f_col)))
})

test_that("spde: Weibull survival with right censoring over an SPDE field", {
  skip_if_not_installed("fmesher")
  set.seed(126)
  n <- 200
  loc <- cbind(runif(n, 0, 1), runif(n, 0, 1))
  mesh <- fmesher::fm_mesh_2d(loc = loc, max.edge = c(0.2, 0.5),
                              cutoff = 0.05)
  fem <- fmesher::fm_fem(mesh)
  fem <- list(c0 = fem$c0, g1 = fem$g1, g2 = fem$g2)
  idx <- mesh$idx$loc

  # Draw the field from the model's own precision, otherwise the SPDE
  # parameters are unidentified and the two optimizers wander apart on a
  # flat ridge instead of disagreeing about the likelihood.
  tau_t <- exp(-1); kap_t <- exp(1.5)
  Q_t <- tau_t^2 * (kap_t^4 * fem$c0 + 2 * kap_t^2 * fem$g1 + fem$g2)
  x_t <- as.vector(backsolve(chol(as.matrix(Q_t)), rnorm(nrow(fem$c0))))

  sexv <- rbinom(n, 1, 0.5)
  agev <- rnorm(n)
  eta_true <- -0.5 + 0.3 * sexv + 0.2 * agev + x_t[idx]
  time <- rweibull(n, shape = 1.2, scale = exp(-eta_true))
  cens <- as.integer(time > 2)          # brms convention: 1 = right censored
  time <- pmin(time, 2)
  dd <- data.frame(time = time, cens = cens, sex = sexv, age = agev,
                   node = factor(idx, levels = seq_len(mesh$n)))

  nll_ref <- function(p) {
    tau <- exp(p$log_tau); kap <- exp(p$log_kappa)
    Q <- tau^2 * (kap^4 * fem$c0 + 2 * kap^2 * fem$g1 + fem$g2)
    om <- exp(p$log_omega)
    eta <- p$beta[1] + p$beta[2] * dd$sex + p$beta[3] * dd$age +
      p$x[idx]
    lam <- exp(eta)
    # ifelse() would evaluate both branches and strip the advector
    # class; index the two groups instead, as the example's loop does
    ob <- which(dd$cens == 0)
    cn <- which(dd$cens == 1)
    ll_ob <- RTMB::dweibull(dd$time[ob], shape = om, scale = 1 / lam[ob],
                            log = TRUE)
    ll_cn <- log(1 - RTMB::pweibull(dd$time[cn], shape = om,
                                    scale = 1 / lam[cn]))
    -RTMB::dgmrf(p$x, 0, Q, log = TRUE) - sum(ll_ob) - sum(ll_cn)
  }
  obj <- RTMB::MakeADFun(nll_ref,
                         list(beta = c(-0.5, 0, 0), log_tau = -2,
                              log_kappa = 1, log_omega = 0,
                              x = rep(0, nrow(fem$c0))),
                         random = "x", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)

  fit <- frm(bf(time | cens(cens) ~ sex + age + spde(fem, gr = node)),
             family = weibull(), data = dd)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-4)
})

test_that("spde() resolves its mesh matrices from data2", {
  # moved from test-data2.R's "every structural special resolves
  # from data2", which keeps the other structural specials
  set.seed(13)
  nn <- 8
  fem <- d2_chain_fem(nn)
  sd_ <- data.frame(node = factor(rep(seq_len(nn), each = 6)),
                    x = stats::rnorm(nn * 6))
  sd_$y <- 0.3 + 0.4 * sd_$x +
    stats::rnorm(nn, 0, 0.7)[as.integer(sd_$node)] +
    stats::rnorm(nrow(sd_), 0, 0.3)
  f_spde <- frm(bf(y ~ x + spde(fem, gr = node)) + gaussian(), data = sd_)
  f_spde2 <- frm(bf(y ~ x + spde(femd2, gr = node)) + gaussian(),
                 data = sd_, data2 = list(femd2 = fem))
  expect_equal(as.numeric(logLik(f_spde2)), as.numeric(logLik(f_spde)),
               tolerance = 1e-10)
})

