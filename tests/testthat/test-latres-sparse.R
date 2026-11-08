# The joint covariance of a latent residual, held sparse (R/joint-sparse.R):
# the readers of the fitted rows and of new rows against the dense route on
# the same fit, the limits that still hold where a dense block is formed, and
# the fallbacks of the row variances.

lrs_data <- function(seed = 11L, G = 30L, K = 12L) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  d$s <- factor(sample.int(6L, n, TRUE))
  e <- as.vector(replicate(G, stats::arima.sim(list(ar = 0.6), K, sd = 0.6)))
  eta <- 0.6 + 0.4 * d$x + e + c(-0.3, 0.2, 0, 0.1, -0.1, 0.3)[d$s]
  d$cnt <- stats::rpois(n, exp(eta))
  d$zc <- ifelse(stats::runif(n) < 0.2, 0L, d$cnt)
  d$ord <- cut(2 * (eta - 0.6) + stats::rlogis(n), c(-Inf, -1, 0, 1, Inf),
               ordered_result = TRUE)
  d
}

# a parameter at its limit is not this file's subject, and a reader of such
# a fit says so in either route
lrs_said <- c("Optimizer did not report convergence",
              "The Laplace approximation of",
              "Large maximum absolute gradient",
              "Standard errors are not available",
              "predictions move along a direction")

lrs_fit <- function(f, family, d) {
  allow_warnings(frm(f, family = family, data = d), lrs_said)
}

# the same fit read through the dense route: its sdreport kept, its joint
# covariance not
lrs_dense <- function(fit) {
  f2 <- fit
  keep <- setdiff(ls(fit$cache), "Vjoint")
  f2$cache <- list2env(mget(keep, envir = fit$cache), parent = emptyenv())
  f2
}

# the dense route, as get_joint_cov() took it before the sparse one
lrs_local_dense <- function(env = parent.frame()) {
  local_mocked_bindings(jcs_joint_cov = function(fit) NULL,
                        .package = "frmtmb", .env = env)
}

lrs_ratio <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  expect_identical(is.finite(a), is.finite(b))
  ok <- is.finite(b)
  if (!any(ok)) return(0)
  max(abs(a[ok] / b[ok] - 1))
}

lrs_count <- function() {
  env <- new.env()
  env$n <- 0L
  env$h <- function(m) {
    env$n <- env$n + 1L
    invokeRestart("muffleMessage")
  }
  env
}

lrs_eps <- sqrt(.Machine$double.eps)

test_that("fitted rows above the old row limit keep their standard errors", {
  d <- lrs_data()
  for (resp in c("cnt", "zc")) local({
    fam <- if (resp == "cnt") poisson() else zero_inflated_poisson()
    fit <- lrs_fit(stats::as.formula(paste(resp, "~ x + (1 | s) +",
                                           "ar(t, g, cov = TRUE)")), fam, d)
    said <- lrs_count()
    got <- withr::with_options(
      list(frmtmb.latres_se_rows = 10L, frmtmb.latres_dense_rows = 10L),
      withCallingHandlers(
        list(v = fitted(fit), r = residuals(fit),
             l = frm_linpred(fit, se.fit = TRUE)),
        frmtmb_latres_se_skipped = said$h))
    expect_identical(said$n, 0L)
    expect_true(all(is.finite(got$v[, "Est.Error"])))
    expect_true(all(is.finite(got$l$se.fit)))
    expect_true(frmtmb:::jc_is_sparse(frmtmb:::get_joint_cov(fit)))
    # the dense route on the same fit, which its default limits let run
    lrs_local_dense()
    fd <- lrs_dense(fit)
    expect_lt(lrs_ratio(got$v[, "Est.Error"], fitted(fd)[, "Est.Error"]),
              lrs_eps)
    expect_lt(lrs_ratio(got$r[, "Est.Error"],
                        residuals(fd)[, "Est.Error"]), lrs_eps)
    expect_lt(lrs_ratio(got$l$se.fit, frm_linpred(fd, se.fit = TRUE)$se.fit),
              lrs_eps)
    expect_false(frmtmb:::jc_is_sparse(frmtmb:::get_joint_cov(fd)))
  })
})

test_that("new rows beside a block keep their SEs above the dense limit", {
  d <- lrs_data(12L)
  fit <- lrs_fit(cnt ~ x + (1 | s) + ar(t, g, cov = TRUE), poisson(), d)
  nd <- data.frame(t = 1:4, g = factor(rep("new", 4L)), x = c(-1, 0, 0.5, 1),
                   s = factor(1:4, levels = levels(d$s)))
  said <- lrs_count()
  got <- withr::with_options(
    list(frmtmb.latres_se_rows = 10L, frmtmb.latres_dense_rows = 10L),
    withCallingHandlers(
      list(v = fitted(fit, newdata = nd),
           l = frm_linpred(fit, newdata = nd, se.fit = TRUE),
           # a curve on new rows reads the estimates and the levels alone
           cu = frm_curve(fit, newdata = nd)),
      frmtmb_latres_se_skipped = said$h))
  expect_identical(said$n, 0L)
  expect_true(all(is.finite(got$v[, "Est.Error"])))
  expect_true(all(is.finite(got$cu$.se)))
  lrs_local_dense()
  fd <- lrs_dense(fit)
  expect_lt(lrs_ratio(got$v[, "Est.Error"],
                      fitted(fd, newdata = nd)[, "Est.Error"]), lrs_eps)
  expect_lt(lrs_ratio(got$l$se.fit,
                      frm_linpred(fd, newdata = nd, se.fit = TRUE)$se.fit),
            lrs_eps)
})

test_that("a fit keeps no covariance of its rows' square after fitted()", {
  d <- lrs_data(13L)
  fit <- lrs_fit(cnt ~ x + (1 | s) + ar(t, g, cov = TRUE), poisson(), d)
  invisible(fitted(fit))
  nb <- length(fit$estimates[["b"]])
  kept <- length(serialize(fit$cache$Vjoint, NULL))
  # the dense joint covariance alone is 8 nb^2 bytes
  expect_lt(kept, 0.25 * 8 * nb^2)
})

test_that("what returns the joint covariance still forms it, under the limit", {
  d <- lrs_data(14L)
  fit <- lrs_fit(cnt ~ x + (1 | s) + ar(t, g, cov = TRUE), poisson(), d)
  jc <- frm_joint_cov(fit)
  expect_true(is.matrix(jc$V))
  lb <- frm_lp_basis(fit)
  expect_true(is.matrix(lb$V))
  expect_identical(dim(lb$V), c(ncol(lb$A), ncol(lb$A)))
  expect_lt(lrs_ratio(sqrt(rowSums((lb$A %*% lb$V) * lb$A)),
                      frm_linpred(fit, se.fit = TRUE)$se.fit), lrs_eps)
  # a new row loads no fitted latent, so none is a column of its design
  nd <- data.frame(t = 1:3, g = factor(rep("new", 3L)), x = 0,
                   s = factor(1:3, levels = levels(d$s)))
  ln <- frm_lp_basis(fit, newdata = nd)
  lat <- frmtmb:::frame_latres(fit$frame)[[1L]][["b_idx"]]
  b_pos <- which(jc$names == "b")
  expect_false(any(ln$coef_pos %in% b_pos[lat]))
  local({
    lrs_local_dense()
    V0 <- frm_joint_cov(lrs_dense(fit))$V
    expect_lt(max(abs(jc$V - V0)) / max(abs(V0)), lrs_eps)
  })
  withr::local_options(frmtmb.latres_dense_rows = 10L)
  f2 <- lrs_dense(fit)
  expect_error(frm_lp_basis(f2), class = "frmtmb_latres_dense")
  expect_error(frm_joint_cov(f2), class = "frmtmb_latres_dense")
  expect_true(all(is.finite(frm_lp_basis(f2, newdata = nd)$V)))
})

test_that("a route with a dense block over the latents skips above the limit", {
  d <- lrs_data(15L)
  fit <- lrs_fit(ord ~ x + ar(t, g, cov = TRUE), cumulative(), d)
  said <- lrs_count()
  quiet <- lrs_count()
  got <- withr::with_options(list(frmtmb.latres_se_rows = 10L), {
    list(v = withCallingHandlers(allow_warnings(fitted(fit), lrs_said),
                                 frmtmb_latres_se_skipped = said$h),
         l = withCallingHandlers(
           allow_warnings(frm_linpred(fit, se.fit = TRUE), lrs_said),
           frmtmb_latres_se_skipped = quiet$h))
  })
  # the category probabilities skip; the predictor, one number per row,
  # reads the factor
  expect_gte(said$n, 1L)
  expect_identical(quiet$n, 0L)
  expect_true(all(is.na(got$v[, "Est.Error", ])))
  expect_true(all(is.finite(got$v[, "Estimate", ])))
  lrs_local_dense()
  ref <- allow_warnings(frm_linpred(lrs_dense(fit), se.fit = TRUE)$se.fit,
                        lrs_said)
  expect_lt(lrs_ratio(got$l$se.fit, ref), lrs_eps)
})

test_that("without the sparse route the skip and the refusal hold", {
  d <- lrs_data(16L)
  fit <- lrs_fit(cnt ~ x + (1 | s) + ar(t, g, cov = TRUE), poisson(), d)
  nd <- data.frame(t = 1:2, g = factor(rep("new", 2L)), x = 0,
                   s = factor(1:2, levels = levels(d$s)))
  lrs_local_dense()
  withr::local_options(frmtmb.latres_se_rows = 10L,
                       frmtmb.latres_dense_rows = 10L)
  said <- lrs_count()
  v <- withCallingHandlers(fitted(lrs_dense(fit)),
                           frmtmb_latres_se_skipped = said$h)
  expect_gte(said$n, 1L)
  expect_true(all(is.na(v[, "Est.Error"])))
  said <- lrs_count()
  v <- withCallingHandlers(fitted(lrs_dense(fit), newdata = nd),
                           frmtmb_latres_se_skipped = said$h)
  expect_gte(said$n, 1L)
  expect_true(all(is.na(v[, "Est.Error"])))
})

test_that("a pair outside the factor's pattern costs one solve per column", {
  # an AR chain's precision, tridiagonal, and five rows that load two latents
  # far apart, whose pair the pattern of the factor does not hold
  n <- 300L
  phi <- 0.6
  Q <- Matrix::bandSparse(n, k = c(0L, 1L), symmetric = TRUE,
                          diagonals = list(c(1, rep(1 + phi^2, n - 2L), 1) +
                                             0.5, rep(-phi, n - 1L)))
  Q <- Matrix::forceSymmetric(methods::as(Q, "CsparseMatrix"))
  S <- list(L = Matrix::Cholesky(Q, perm = TRUE, super = TRUE, LDL = FALSE),
            r = seq_len(n))
  Ar <- Matrix::sparseMatrix(i = c(seq_len(n), 1:5),
                             j = c(seq_len(n), n - 0:4),
                             x = c(rep(1, n), rep(0.5, 5L)), dims = c(n, n))
  Ad <- as.matrix(Ar)
  ref <- rowSums((Ad %*% solve(as.matrix(Q))) * Ad)
  # the five far pairs are not in the pattern, so a lookup alone misses them
  keys <- frmtmb:::jcs_selinv_keys(S)
  far <- match((0:4) * keys$N + (n - 1L - 0:4), keys$key)
  expect_true(all(is.na(far)))
  cols <- integer(0)
  pairs <- frmtmb:::jcs_W_pairs
  got <- local({
    local_mocked_bindings(
      jcs_W_pairs = function(S, i, j) {
        v <- pairs(S, i, j)
        cols <<- c(cols, attr(v, "columns"))
        v
      },
      jcs_rowvar_solve = function(S, Ar, ri) stop("one solve per row"),
      .package = "frmtmb")
    frmtmb:::jcs_rowvar_inner(S, Ar, seq_len(n))
  })
  expect_lt(max(abs(got / ref - 1)), lrs_eps)
  # five pairs, five columns at most, and no row solved on its own
  expect_identical(length(cols), 1L)
  expect_lte(cols, 5L)
  # where TMB does not give the selected inverse, every row is solved
  local_mocked_bindings(jcs_selinv = function(S) NULL, .package = "frmtmb")
  expect_lt(max(abs(frmtmb:::jcs_rowvar_inner(S, Ar, seq_len(n)) / ref - 1)),
            lrs_eps)
})

# A hurdle or zero-inflated mean reads mu's latent and a level of `hu` or `zi`
# on one row, a pair no row's likelihood couples, so it is outside the
# factor's pattern on every row (hurdle) or on each row with y > 0 (zi).
lrs_split_data <- function(seed, n, hurdle) {
  set.seed(seed)
  d <- data.frame(t = seq_len(n), g = factor(1), x = stats::rnorm(n),
                  id = factor(sample.int(12L, n, TRUE)))
  e <- as.numeric(stats::arima.sim(list(ar = 0.5), n, sd = 0.5))
  u <- stats::rnorm(12L, 0, 0.5)
  mu <- exp(0.5 + 0.3 * d$x + e)
  y <- stats::rpois(n, mu)
  if (hurdle) {
    while (any(y == 0)) y[y == 0] <- stats::rpois(sum(y == 0), mu[y == 0])
  }
  d$y <- ifelse(stats::runif(n) < stats::plogis(-1 + u[d$id]), 0L, y)
  d
}

test_that("a hurdle's or a zi's pairs cost one solve per level, not per row", {
  for (hurdle in c(TRUE, FALSE)) local({
    d <- lrs_split_data(21L, 400L, hurdle)
    fit <- if (hurdle) {
      lrs_fit(bf(y ~ x + ar(t, g, cov = TRUE), hu ~ (1 | id)),
              hurdle_poisson(), d)
    } else {
      lrs_fit(bf(y ~ x + ar(t, g, cov = TRUE), zi ~ (1 | id)),
              zero_inflated_poisson(), d)
    }
    cols <- integer(0)
    pairs <- frmtmb:::jcs_W_pairs
    v <- local({
      local_mocked_bindings(
        jcs_W_pairs = function(S, i, j) {
          out <- pairs(S, i, j)
          cols <<- c(cols, attr(out, "columns"))
          out
        },
        jcs_rowvar_solve = function(S, Ar, ri) stop("one solve per row"),
        .package = "frmtmb")
      allow_warnings(fitted(fit), lrs_said)
    })
    # the pairs left the pattern, and took one column per level of the gate
    expect_gte(length(cols), 1L)
    expect_lte(max(cols), nlevels(d$id))
    lrs_local_dense()
    ref <- allow_warnings(fitted(lrs_dense(fit)), lrs_said)
    expect_lt(lrs_ratio(v[, "Est.Error"], ref[, "Est.Error"]), lrs_eps)
  })
})

test_that("the standard-error check's propagation covariance is the one read", {
  # no group effect in the data, so (1 | id) runs to its boundary, its sd
  # loses its standard error, and the check keeps a covariance for
  # propagation (the repair branch); the reference forms joint_cov_repair()'s
  # block formula from it here (helper-latres-ref.R)
  set.seed(20261106L)
  d <- expand.grid(t = seq_len(14L), g = factor(seq_len(30L)))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  d$id <- factor(sample.int(8L, n, TRUE))
  e <- as.vector(replicate(30L, stats::arima.sim(list(ar = 0.6), 14L,
                                                 sd = 0.5)))
  d$y <- stats::rpois(n, exp(0.5 + 0.4 * d$x + e))
  fit <- suppressMessages(lrs_fit(y ~ x + (1 | id) + ar(t, g, cov = TRUE),
                                  poisson(), d))
  # the point is the repair branch: a fit that stops inside the boundary
  # keeps no propagation covariance and would test solve(Q) instead, so it
  # fails here rather than skip
  expect_false(is.null(sdr_of(fit)$cov_fixed_prop))
  ref <- latres_ref_se(fit)
  got <- suppressMessages(frm_linpred(fit, se.fit = TRUE))$se.fit
  # agreement is to rounding: 4.4e-16 here, where reading solve(Q) in place
  # of the kept covariance moves it by 2.7e-10
  tol <- 1000 * .Machine$double.eps
  expect_lt(latres_ref_ratio(got, ref), tol)
  V <- suppressMessages(frm_joint_cov(fit))$V
  Vr <- latres_ref_V(fit)
  ok <- is.finite(V)
  expect_lt(max(abs(V[ok] - Vr[ok])) / max(abs(Vr)), tol)
})
