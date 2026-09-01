# Heavy reference-validation file: full runs happen locally and in CI
# (NOT_CRAN=true). CRAN only needs to see the package work, not the
# validation program re-derived, and this file is a main contributor
# to CRAN-condition check time.
skip_on_cran()

# v0.12: function-on-function and af-style functional terms, sparse
# precision gr(), new-level prediction variance.

test_that("function-on-function regression matches mgcv exactly", {
  set.seed(801)
  n_sub <- 60; Ks <- 20; Kt <- 15
  ss <- seq(0, 1, length.out = Ks)
  tt <- seq(0, 1, length.out = Kt)
  Xf <- t(replicate(n_sub, cumsum(rnorm(Ks, 0, 0.4))))
  B <- outer(ss, tt, function(s, t) 2 * sin(pi * s) * cos(pi * t))
  Ymat <- Xf %*% B / Ks + matrix(rnorm(n_sub * Kt, 0, 0.2), n_sub)

  dd <- expand.grid(t_idx = seq_len(Kt), i = seq_len(n_sub))
  dd$y <- as.vector(t(Ymat))
  n_row <- nrow(dd)
  dd$Smat <- matrix(ss, n_row, Ks, byrow = TRUE)
  dd$Tmat <- matrix(tt[dd$t_idx], n_row, Ks)
  dd$Lmat <- Xf[dd$i, ] / Ks

  fit <- frm(bf(y ~ t2(Smat, Tmat, by = Lmat)) + gaussian(), data = dd)
  ref <- mgcv::gam(y ~ t2(Smat, Tmat, by = Lmat), data = dd,
                   method = "ML")
  expect_lt(abs(as.numeric(logLik(fit)) - (-as.numeric(ref$gcv.ubre))),
            1e-3)
  expect_lt(max(abs(fitted(fit) - fitted(ref))), 1e-4)
})

test_that("af-style nonlinear functional terms fit (surfaces match mgcv)", {
  set.seed(802)
  n <- 300; K <- 25
  tg <- seq(0, 1, length.out = K)
  Xc <- t(replicate(n, cumsum(rnorm(K, 0, 0.3))))
  y <- 1 + rowSums(0.5 * Xc^2 - Xc * matrix(tg, n, K, byrow = TRUE)) / K +
    rnorm(n, 0, 0.3)
  dd <- data.frame(y = y)
  dd$Tmat <- matrix(tg, n, K, byrow = TRUE)
  dd$Xmat <- Xc
  dd$Wmat <- matrix(1 / K, n, K)

  fit <- suppressWarnings(
    frm(bf(y ~ t2(Tmat, Xmat, by = Wmat)) + gaussian(), data = dd)
  )
  ref <- mgcv::gam(y ~ t2(Tmat, Xmat, by = Wmat), data = dd,
                   method = "ML")
  # multi-penalty tensor smoothing profiles are mildly multimodal: the
  # two optimizers land on nearby local modes; the fitted surfaces agree
  expect_lt(abs(as.numeric(logLik(fit)) - (-as.numeric(ref$gcv.ubre))),
            0.5)
  expect_lt(max(abs(fitted(fit) - fitted(ref))), 0.1)
})

test_that("gr(prec=) matches gr(cov=) with the inverse matrix", {
  set.seed(811)
  ng <- 40
  # tri-diagonal (AR1-type) sparse precision
  Q <- Matrix::bandSparse(ng, k = c(-1, 0, 1),
                          diagonals = list(rep(-0.4, ng - 1),
                                           rep(1.2, ng),
                                           rep(-0.4, ng - 1)))
  dimnames(Q) <- list(as.character(seq_len(ng)),
                      as.character(seq_len(ng)))
  A <- solve(as.matrix(Q))
  dimnames(A) <- dimnames(Q)

  b_true <- drop(crossprod(chol(A), rnorm(ng))) * 0.7
  n <- 400
  g <- factor(rep(seq_len(ng), each = n / ng))
  x <- rnorm(n)
  dd <- data.frame(y = 1 + 0.5 * x + b_true[g] + rnorm(n, 0, 0.5),
                   x = x, g = g)

  fp <- frm(bf(y ~ x + (1 | gr(g, prec = Q))) + gaussian(), data = dd)
  fc <- frm(bf(y ~ x + (1 | gr(g, cov = A))) + gaussian(), data = dd)
  expect_lt(abs(as.numeric(logLik(fp)) - as.numeric(logLik(fc))), 1e-6)
  expect_vector_equal(fixef(fp)$mu, fixef(fc)$mu, tol = 1e-5)
  # marginal draws respect the precision structure
  set.seed(1)
  b1 <- frmtmb:::draw_b(fp)
  expect_length(b1, ng)
})

test_that("gr(prec=) takes correlated slopes", {
  # inv(A (x) Sigma) = A^-1 (x) Sigma^-1, so the precision side of the
  # Kronecker product has to reproduce the dense gr(cov=) fit exactly.
  # The design needs real slope variation: a data set whose slope
  # variance is zero puts both fits on the boundary, where they still
  # agree but nothing about the Kronecker assembly is exercised.
  set.seed(3)
  ng <- 8
  Q <- Matrix::bandSparse(ng, k = c(-1, 0, 1),
                          diagonals = list(rep(-0.7, ng - 1),
                                           rep(2, ng),
                                           rep(-0.7, ng - 1)))
  dimnames(Q) <- list(paste0("g", seq_len(ng)), paste0("g", seq_len(ng)))
  A <- solve(as.matrix(Q))
  dimnames(A) <- dimnames(Q)
  S <- matrix(c(0.8, 0.3, 0.3, 0.5), 2, 2)
  b_true <- drop(crossprod(chol(kronecker(A, S)), rnorm(2 * ng)))
  dd <- data.frame(g = factor(rep(rownames(A), each = 12),
                              levels = rownames(A)))
  dd$x <- rnorm(nrow(dd))
  j <- as.integer(dd$g)
  dd$y <- 0.5 + b_true[(j - 1) * 2 + 1] + b_true[(j - 1) * 2 + 2] * dd$x +
    rnorm(nrow(dd), 0, 0.4)

  fp <- frm(bf(y ~ x + (1 + x | gr(g, prec = Q))) + gaussian(), data = dd)
  fc <- frm(bf(y ~ x + (1 + x | gr(g, cov = A))) + gaussian(), data = dd)
  expect_lt(abs(as.numeric(logLik(fp)) - as.numeric(logLik(fc))), 1e-8)
  expect_vector_equal(fp$estimates$theta, fc$estimates$theta, tol = 1e-8)
  expect_vector_equal(fixef(fp)$mu, fixef(fc)$mu, tol = 1e-8)
  expect_vector_equal(VarCorr(fp)[[1]], VarCorr(fc)[[1]], tol = 1e-8)
  expect_vector_equal(ranef(fp)[[1]], ranef(fc)[[1]], tol = 1e-8)
  nd <- dd[c(1L, 30L, 90L), ]
  expect_vector_equal(predict(fp, newdata = nd, se.fit = TRUE)$se.fit,
                      predict(fc, newdata = nd, se.fit = TRUE)$se.fit,
                      tol = 1e-8)
  set.seed(1)
  expect_length(frmtmb:::draw_b(fp), 2 * ng)
})

test_that("new levels add the block variance to prediction SEs", {
  set.seed(812)
  n <- 300
  dd <- data.frame(x = rnorm(n), g = factor(rep(1:15, 20)))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(15, 0, 0.8)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

  nd <- data.frame(x = 0, g = factor("NEW"))
  p_new <- predict(fit, newdata = nd, se.fit = TRUE,
                   allow_new_levels = TRUE)
  p_pop <- predict(fit, newdata = data.frame(x = 0), re.form = NA,
                   se.fit = TRUE)
  # same point prediction, inflated uncertainty
  expect_equal(p_new$fit, p_pop$fit, tolerance = 1e-8)
  sd_g2 <- VarCorr(fit)[[1]][1, 1]
  expect_equal(p_new$se.fit^2, p_pop$se.fit^2 + sd_g2, tolerance = 1e-6)

  # known levels are unaffected
  nd2 <- data.frame(x = 0, g = factor("3", levels = levels(dd$g)))
  p_known <- predict(fit, newdata = nd2, se.fit = TRUE)
  expect_lt(p_known$se.fit, p_new$se.fit)
})