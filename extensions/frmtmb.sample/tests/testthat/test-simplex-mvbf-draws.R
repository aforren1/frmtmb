# Draws of a logistic_normal() response inside mvbf(): the correlation
# is reported as lncor_<resp>__b__c, as core's variables() names it, and
# it takes brms's lkj(1) default under the response's name. frmtmb
# 0.71.0 refused the model.

test_that("a multivariate logistic_normal() samples under its own names", {
  skip_sampler()
  set.seed(7)
  n <- 200
  dd <- data.frame(x = stats::rnorm(n))
  W <- matrix(stats::rnorm(2 * n), n) %*% chol(matrix(c(1, .4, .4, 1), 2))
  E <- cbind(0, 0.4 + 0.5 * dd$x + 0.6 * W[, 1],
             -0.3 - 0.4 * dd$x + 0.9 * W[, 2])
  Y <- exp(E) / rowSums(exp(E))
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd$w <- 1 + 0.5 * dd$x + stats::rnorm(n)
  fit <- frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                  bf(w ~ x, family = gaussian())), data = dd)
  msg <- character()
  allow_warnings(
    ds <- withCallingHandlers(
      frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 2),
      message = function(m) {
        msg <<- c(msg, conditionMessage(m))
        invokeRestart("muffleMessage")
      }),
    c("ESS", "R-hat"))
  expect_true(any(grepl("lncor (resp = Y)   lkj(1)", msg, fixed = TRUE)))
  v <- variables(ds)
  expect_identical(grep("lncor", v, value = TRUE), "lncor_Y__b__c")
  # every draw is a correlation, and the draws center near the estimate
  r <- as.matrix(posterior_summary(ds))["lncor_Y__b__c", ]
  est <- summary(fit)$spec_pars["lncor_Y(b,c)", ]
  expect_true(r[["Q2.5"]] > -1 && r[["Q97.5"]] < 1)
  expect_lt(abs(r[["Estimate"]] - est[["Estimate"]]),
            4 * est[["Est.Error"]])
})

test_that("the draws name the correlation by brms's Stan name", {
  # Y_a is brms's Ya (review finding 1); the first build wrote lncor_Y_a
  skip_sampler()
  set.seed(8)
  n <- 200
  dd <- data.frame(x = stats::rnorm(n))
  G <- matrix(stats::rgamma(3 * n, shape = 4), n, 3)
  Y <- G / rowSums(G)
  colnames(Y) <- c("a", "b", "c")
  dd$Y_a <- Y
  dd$w <- stats::rnorm(n)
  fit <- frm(mvbf(bf(Y_a ~ x, family = logistic_normal()),
                  bf(w ~ x, family = gaussian())), data = dd)
  allow_warnings(
    ds <- suppressMessages(frm_sample(fit, chains = 1, iter = 200,
                                      refresh = 0, seed = 2)),
    c("ESS", "R-hat"))
  expect_identical(grep("lncor", variables(ds), value = TRUE),
                   "lncor_Ya__b__c")
})
