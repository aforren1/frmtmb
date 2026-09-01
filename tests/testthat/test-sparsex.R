# frmtmb_control(sparse_x = TRUE): parametric fixed-effect designs as
# sparse Matrix objects for many-level fixed factors (glmmTMB sparseX).

test_that("sparse_x reproduces the dense fit for a many-level factor", {
  set.seed(7)
  n <- 4000
  dd <- data.frame(
    f = factor(rep(sprintf("s%03d", 1:200), each = 20)),
    x = rnorm(n), z = rnorm(n),
    g = factor(sample(1:25, n, TRUE))
  )
  fe <- rnorm(200, 0, 0.8)
  dd$y <- rnorm(n, 1 + fe[as.integer(dd$f)] + 0.5 * dd$x - 0.3 * dd$z +
                  rnorm(25, 0, 0.4)[dd$g], 1)
  form <- bf(y ~ f + x + z + (1 | g)) + gaussian()
  ctl <- function(sparse) {
    frmtmb_control(sparse_x = sparse, restarts = 2,
                   optCtrl = list(iter.max = 2000, eval.max = 2000,
                                  rel.tol = 1e-12))
  }
  f_d <- frm(form, data = dd, control = ctl(FALSE))
  f_s <- frm(form, data = dd, control = ctl(TRUE))

  X_d <- f_d$frame$linpreds[["y.mu"]]$X
  X_s <- f_s$frame$linpreds[["y.mu"]]$X
  expect_s4_class(X_s, "dgCMatrix")
  expect_identical(colnames(X_s), colnames(X_d))
  # memory: the sparse design is a small fraction of the dense one
  expect_lt(as.numeric(object.size(X_s)) / as.numeric(object.size(X_d)),
            0.25)

  # the sparse and dense tapes are the same function
  expect_equal(as.numeric(f_s$obj$fn(f_d$opt$par)),
               as.numeric(f_d$opt$objective), tolerance = 1e-10)
  expect_loglik_equal(f_s, f_d, tol = 1e-8)
  # independent optimizer runs stop within trajectory noise of each other
  expect_vector_equal(unlist(fixef(f_s)), unlist(fixef(f_d)), tol = 1e-5)

  # polishing the sparse fit from the dense optimum lands on the
  # identical solution: the sparse objective has the same optimum
  st <- list(beta = unname(f_d$estimates$beta),
             betad = unname(f_d$estimates$betad),
             theta = unname(f_d$estimates$theta))
  f_w <- suppressWarnings(frm(form, data = dd, control = ctl(TRUE),
                              start = st))
  expect_loglik_equal(f_w, f_d, tol = 1e-8)
  expect_vector_equal(unlist(fixef(f_w)), unlist(fixef(f_d)), tol = 1e-8)
  expect_vector_equal(diag(vcov(f_w)), diag(vcov(f_d)), tol = 1e-8)

  # method surface: predict/se.fit (in-sample and newdata), summary,
  # confint, model.matrix
  p_d <- predict(f_d, se.fit = TRUE)
  p_w <- predict(f_w, se.fit = TRUE)
  expect_vector_equal(p_w$fit, p_d$fit, tol = 1e-8)
  expect_vector_equal(p_w$se.fit, p_d$se.fit, tol = 1e-8)
  nd <- dd[seq_len(40), ]
  pn_d <- predict(f_d, newdata = nd, se.fit = TRUE)
  pn_w <- predict(f_w, newdata = nd, se.fit = TRUE)
  expect_vector_equal(pn_w$fit, pn_d$fit, tol = 1e-8)
  expect_vector_equal(pn_w$se.fit, pn_d$se.fit, tol = 1e-8)
  expect_identical(rownames(summary(f_w)$coefficients$mu),
                   rownames(summary(f_d)$coefficients$mu))
  ci_d <- confint(f_d, parm = c("x", "z"))
  ci_w <- confint(f_w, parm = c("x", "z"))
  expect_vector_equal(ci_w, ci_d, tol = 1e-8)
  expect_s4_class(model.matrix(f_w), "dgCMatrix")
})

test_that("sparse_x is exact for a distributional (sigma ~ x) model", {
  set.seed(8)
  dd <- data.frame(x = rnorm(400), w = rnorm(400),
                   g = factor(rep(1:20, each = 20)))
  dd$y <- rnorm(400, 1 + 0.6 * dd$x + rnorm(20, 0, 0.5)[dd$g],
                exp(-0.2 + 0.3 * dd$x))
  form <- bf(y ~ x + w + (1 | g), sigma ~ x) + gaussian()
  ctl <- function(sparse) frmtmb_control(sparse_x = sparse, restarts = 2)
  f_d <- frm(form, data = dd, control = ctl(FALSE))
  f_s <- frm(form, data = dd, control = ctl(TRUE))
  expect_s4_class(f_s$frame$linpreds[["y.sigma"]]$X, "dgCMatrix")
  expect_loglik_equal(f_s, f_d, tol = 1e-8)
  expect_vector_equal(unlist(fixef(f_s)), unlist(fixef(f_d)), tol = 1e-8)
  expect_vector_equal(diag(vcov(f_s)), diag(vcov(f_d)), tol = 1e-8)
  ps_d <- predict(f_d, dpar = "sigma", se.fit = TRUE)
  ps_s <- predict(f_s, dpar = "sigma", se.fit = TRUE)
  expect_vector_equal(ps_s$fit, ps_d$fit, tol = 1e-8)
  expect_vector_equal(ps_s$se.fit, ps_d$se.fit, tol = 1e-8)
})

test_that("sparse_x is exact for mo() terms (patched zero columns)", {
  set.seed(9)
  dd <- data.frame(x = rnorm(500), m = sample(0:4, 500, TRUE),
                   g = factor(rep(1:25, each = 20)))
  dd$y <- rnorm(500, 1 + 0.5 * dd$x +
                  2 * cumsum(c(0, 0.4, 0.3, 0.2, 0.1))[dd$m + 1L] +
                  rnorm(25, 0, 0.3)[dd$g], 1)
  form <- bf(y ~ x + mo(m) + (1 | g)) + gaussian()
  f_d <- suppressWarnings(frm(form, data = dd))
  f_s <- suppressWarnings(
    frm(form, data = dd, control = frmtmb_control(sparse_x = TRUE)))
  expect_loglik_equal(f_s, f_d, tol = 1e-8)
  expect_vector_equal(unlist(fixef(f_s)), unlist(fixef(f_d)), tol = 1e-8)
  p_d <- predict(f_d, se.fit = TRUE)
  p_s <- predict(f_s, se.fit = TRUE)
  expect_vector_equal(p_s$fit, p_d$fit, tol = 1e-8)
  expect_vector_equal(p_s$se.fit, p_d$se.fit, tol = 1e-8)
  nd <- dd[seq_len(20), ]
  pn_d <- predict(f_d, newdata = nd, se.fit = TRUE)
  pn_s <- predict(f_s, newdata = nd, se.fit = TRUE)
  expect_vector_equal(pn_s$fit, pn_d$fit, tol = 1e-8)
  expect_vector_equal(pn_s$se.fit, pn_d$se.fit, tol = 1e-8)
})

test_that("sparse_x drops the same rank-deficient columns as dense", {
  set.seed(10)
  dd <- data.frame(x = rnorm(200), g = factor(rep(1:10, 20)))
  dd$x2 <- 2 * dd$x
  dd$y <- rnorm(200, 1 + dd$x, 1)
  form <- bf(y ~ x + x2 + (1 | g)) + gaussian()
  expect_message(f_d <- frm(form, data = dd), "rank deficient")
  expect_message(
    f_s <- frm(form, data = dd,
               control = frmtmb_control(sparse_x = TRUE)),
    "dropping column\\(s\\): x2")
  expect_identical(colnames(f_s$frame$linpreds[["y.mu"]]$X),
                   colnames(f_d$frame$linpreds[["y.mu"]]$X))
  expect_vector_equal(unlist(fixef(f_s)), unlist(fixef(f_d)), tol = 1e-8)
})

test_that("sparse_x keeps dense naming, smooths, and NA semantics", {
  set.seed(11)
  n <- 400
  dd <- data.frame(x = runif(n, -2, 2), g = factor(rep(1:20, each = 20)),
                   f = factor(sample(letters[1:8], n, TRUE)))
  dd$y <- rnorm(n, 1 + sin(2 * dd$x) + 0.3 * as.integer(dd$f) +
                  rnorm(20, 0, 0.4)[dd$g], 0.5)

  # matrix-valued terms keep the dense poly(x, 2)k names
  f_p <- frm(bf(y ~ poly(x, 2) + (1 | g)) + gaussian(), data = dd,
             control = frmtmb_control(sparse_x = TRUE))
  expect_identical(names(fixef(f_p)$mu),
                   c("(Intercept)", "poly(x, 2)1", "poly(x, 2)2"))

  # smooth null-space columns cbind onto the sparse X; ML and REML
  form <- bf(y ~ f + s(x) + (1 | g)) + gaussian()
  f_d <- frm(form, data = dd)
  f_s <- frm(form, data = dd, control = frmtmb_control(sparse_x = TRUE))
  expect_loglik_equal(f_s, f_d, tol = 1e-8)
  expect_vector_equal(unlist(fixef(f_s)), unlist(fixef(f_d)), tol = 1e-8)
  nd <- data.frame(x = seq(-2, 2, length.out = 25),
                   f = factor("c", levels = letters[1:8]),
                   g = factor(1, levels = levels(dd$g)))
  p_d <- predict(f_d, newdata = nd, se.fit = TRUE, re.form = NA)
  p_s <- predict(f_s, newdata = nd, se.fit = TRUE, re.form = NA)
  expect_vector_equal(p_s$fit, p_d$fit, tol = 1e-8)
  expect_vector_equal(p_s$se.fit, p_d$se.fit, tol = 1e-8)
  r_d <- frm(form, data = dd, REML = TRUE)
  r_s <- frm(form, data = dd, REML = TRUE,
             control = frmtmb_control(sparse_x = TRUE))
  expect_loglik_equal(r_s, r_d, tol = 1e-8)
  expect_vector_equal(unlist(fixef(r_s)), unlist(fixef(r_d)), tol = 1e-8)

  # NA factor levels in newdata propagate to NA predictions (the sparse
  # builder zeroes NA factor rows, so those frames fall back to dense)
  ndna <- data.frame(x = c(0, 0.5), f = factor(c("a", NA),
                                               levels = letters[1:8]),
                     g = factor(c(1, 2), levels = levels(dd$g)))
  pna_d <- predict(f_d, newdata = ndna, re.form = NA)
  pna_s <- predict(f_s, newdata = ndna, re.form = NA)
  expect_identical(is.na(pna_s), is.na(pna_d))
  expect_true(is.na(pna_s[2]))
  expect_equal(pna_s[1], pna_d[1], tolerance = 1e-8)
})
