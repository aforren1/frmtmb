# The draws surface on core's simplex families, dirichlet() and
# logistic_normal(), in brms 2.23.0's shapes (read off brms fits in
# dev/famsimplex-brms-fits2.R, recorded in dev/famsimplex-findings.md).
# log_lik() is checked against densities written here from lgamma(),
# log() and a hand multivariate normal, not from the package code. Each
# reference carries the sum of its absolute terms as attribute "scale",
# which is what rounding is measured against: the lgamma() terms cancel.

sxd_ref_dirichlet <- function(Y, A) {
  t1 <- lgamma(rowSums(A))
  t2 <- lgamma(A)
  t3 <- (A - 1) * log(Y)
  structure(t1 - rowSums(t2) + rowSums(t3),
            scale = abs(t1) + rowSums(abs(t2)) + rowSums(abs(t3)))
}

sxd_ref_ln <- function(Y, M, S, r) {
  X <- log(Y[, -1L, drop = FALSE]) - log(Y[, 1L])
  Z <- (X - M) / S
  q <- (Z[, 1]^2 - 2 * r * Z[, 1] * Z[, 2] + Z[, 2]^2) / (1 - r^2)
  structure(-log(2 * pi) - 0.5 * log(1 - r^2) - 0.5 * q - rowSums(log(S)) -
              rowSums(log(Y)),
            scale = log(2 * pi) + 0.5 * abs(log(1 - r^2)) + 0.5 * q +
              rowSums(abs(log(S))) + rowSums(abs(log(Y))))
}

sxd_data <- function(seed, n = 150) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  E <- cbind(0, 0.4 + 0.6 * d$x, -0.3 - 0.5 * d$x)
  P <- exp(E) / rowSums(exp(E))
  G <- matrix(stats::rgamma(3 * n, shape = 8 * P), n)
  Yd <- G / rowSums(G)
  W <- matrix(stats::rnorm(2 * n), n) %*% chol(matrix(c(1, .3, .3, 1), 2))
  E2 <- cbind(0, E[, 2] + 0.6 * W[, 1], E[, 3] + 0.8 * W[, 2])
  Yl <- exp(E2) / rowSums(exp(E2))
  colnames(Yd) <- colnames(Yl) <- c("a", "b", "c")
  d$Yd <- Yd
  d$Yl <- Yl
  d
}

sxd_draws <- local({
  cache <- list()
  function(which) {
    skip_sampler()
    if (is.null(cache[[which]])) {
      d <- sxd_data(71)
      fit <- switch(which,
        dirichlet = frm(bf(Yd ~ x), family = dirichlet(), data = d),
        ln = frm(bf(Yl ~ x), family = logistic_normal(), data = d))
      msg <- character()
      ds <- withCallingHandlers(
        suppressWarnings(frm_sample(fit, chains = 1, iter = 300,
                                    refresh = 0, seed = 3)),
        message = function(m) {
          msg <<- c(msg, conditionMessage(m))
          invokeRestart("muffleMessage")
        })
      cache[[which]] <<- list(d = d, fit = fit, ds = ds, msg = msg)
    }
    cache[[which]]
  }
})

test_that("log_lik() is the reference density at the draws", {
  for (which in c("dirichlet", "ln")) {
    cs <- sxd_draws(which)
    ll <- log_lik(cs$ds)
    expect_identical(dim(ll), c(150L, 150L))
    idx <- frmtmb.sample:::draws_par_index(cs$ds$fit)
    for (s in c(1L, nrow(ll))) {
      f <- frmtmb.sample:::draws_fit_at(cs$ds, s, idx)
      dp <- frmtmb:::eval_dpars(f)[[1]]
      ref <- if (identical(which, "dirichlet")) {
        E <- cbind(0, dp$mub, dp$muc)
        sxd_ref_dirichlet(cs$d$Yd, exp(E) / rowSums(exp(E)) * dp$phi)
      } else {
        r <- frmtmb:::lncor_matrix(f$estimates$lncor_raw, 2L)[1, 2]
        sxd_ref_ln(cs$d$Yl, cbind(dp$mub, dp$muc),
                   cbind(rep(dp$sigmab, length.out = 150),
                         rep(dp$sigmac, length.out = 150)), r)
      }
      expect_lt(max(abs(ll[s, ] - as.numeric(ref)) / attr(ref, "scale")),
                64 * .Machine$double.eps, label = which)
    }
  }
})

test_that("posterior_epred() and posterior_predict() take brms's shapes", {
  cs <- sxd_draws("dirichlet")
  ep <- posterior_epred(cs$ds, ndraws = 5)
  # brms 2.23.0: 5 x 300 x 3, dimnames list(NULL, NULL, c("a", "b", "c"))
  expect_identical(dim(ep), c(5L, 150L, 3L))
  expect_identical(dimnames(ep)[[3]], c("a", "b", "c"))
  expect_lt(max(abs(apply(ep, 1:2, sum) - 1)), 64 * .Machine$double.eps)
  for (which in c("dirichlet", "ln")) {
    pp <- posterior_predict(sxd_draws(which)$ds, ndraws = 5)
    expect_identical(dim(pp), c(5L, 150L, 3L), label = which)
    expect_identical(dimnames(pp)[[3]], c("a", "b", "c"), label = which)
    expect_lt(max(abs(apply(pp, 1:2, sum) - 1)), 64 * .Machine$double.eps,
              label = which)
  }
  # brms refuses the expectation of a logistic-normal composition
  expect_error(posterior_epred(sxd_draws("ln")$ds, ndraws = 5),
               "Cannot compute expected values of the posterior predictive")
})

test_that("the draws carry brms's lncor__b__c, read back exactly", {
  cs <- sxd_draws("ln")
  v <- variables(cs$ds)
  expect_true("lncor__b__c" %in% v)
  expect_false("lncor_raw_1" %in% v)
  r <- as_draws_matrix(cs$ds, variable = "lncor__b__c")
  expect_true(all(abs(r) < 1))
  # the column is the correlation of the draw the model is evaluated at
  idx <- frmtmb.sample:::draws_par_index(cs$ds$fit)
  f <- frmtmb.sample:::draws_fit_at(cs$ds, 1L, idx)
  r1 <- frmtmb:::lncor_matrix(f$estimates$lncor_raw, 2L)[1, 2]
  expect_lt(abs(r1 - r[1, 1]), 64 * .Machine$double.eps)
  # and the defaults are brms's: lkj(1) on it, student_t on sigma<k>
  expect_true(any(grepl("lncor", cs$msg) & grepl("lkj(1)", cs$msg,
                                                  fixed = TRUE)))
  expect_true(any(grepl("sigmab", cs$msg) & grepl("student_t(3, 0, 2.5)",
                                                   cs$msg, fixed = TRUE)))
  # a dirichlet fit names its phi among the defaults it cannot carry
  expect_true(any(grepl("no defaults for phi", sxd_draws("dirichlet")$msg,
                        fixed = TRUE)))
})

test_that("pp_check() and predictive_error() are refused as in brms", {
  cs <- sxd_draws("dirichlet")
  skip_if_not_installed("bayesplot")
  expect_error(pp_check(cs$ds), "'pp_check' is not implemented")
  expect_error(predictive_error(cs$ds), "not defined for the 'dirichlet'")
})
