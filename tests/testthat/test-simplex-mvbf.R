# The simplex families inside mvbf(). Without a residual correlation the
# multivariate log-likelihood is the sum of the univariate ones, which
# is the reference here; test-refcat-brms.R checks both models against
# brms 2.23.0's Stan program. 0.71.0 refused logistic_normal() in a
# multivariate model ("Families with extra parameters").

smv_data <- function(seed, n = 300) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  W <- matrix(stats::rnorm(2 * n), n) %*% chol(matrix(c(1, .4, .4, 1), 2))
  E <- cbind(0, 0.4 + 0.5 * dd$x + 0.6 * W[, 1],
             -0.3 - 0.4 * dd$x + 0.9 * W[, 2])
  Y <- exp(E) / rowSums(exp(E))
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  G <- matrix(stats::rgamma(n * 3, shape = 4), n, 3)
  V <- G / rowSums(G)
  colnames(V) <- c("u", "v", "w")
  dd$V <- V
  dd$w <- 1 + 0.5 * dd$x + stats::rnorm(n)
  dd
}

smv_ctl <- frmtmb_control(grad_tol = 1e-6, restarts = 3)

# the largest difference, relative to the larger magnitude
rc_rel_smv <- function(a, b) max(abs(a - b)) / max(abs(a), abs(b))

test_that("logistic_normal() in mvbf() is the sum of its two models", {
  dd <- smv_data(1)
  fm <- frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                 bf(w ~ x, family = gaussian())), data = dd,
            control = smv_ctl)
  f1 <- frm(bf(Y ~ x), family = logistic_normal(), data = dd,
            control = smv_ctl)
  f2 <- frm(bf(w ~ x), data = dd, control = smv_ctl)
  ll <- as.numeric(logLik(f1)) + as.numeric(logLik(f2))
  expect_lt(abs(as.numeric(logLik(fm)) - ll), 1e-8 * abs(ll))
  # the correlation is reported per response, lncor_<resp>, beside the
  # univariate fit's lncor(b,c)
  sp <- summary(fm)$spec_pars
  expect_identical(rownames(sp),
                   c("sigmab_Y", "sigmac_Y", "sigma_w", "lncor_Y(b,c)"))
  r1 <- summary(f1)$spec_pars["lncor(b,c)", ]
  rm <- sp["lncor_Y(b,c)", ]
  expect_lt(abs(rm[["Estimate"]] - r1[["Estimate"]]),
            1e-3 * r1[["Est.Error"]])
  expect_lt(abs(rm[["Est.Error"]] / r1[["Est.Error"]] - 1), 1e-3)
  expect_true("lncor_Y__b__c" %in% variables(fm))
  h <- hypothesis(fm, "lncor_Y__b__c = 0", class = NULL)
  expect_lt(abs(h$hypothesis$Estimate - r1[["Estimate"]]),
            1e-3 * r1[["Est.Error"]])
  # brms's class lncor, with the response named
  dp <- default_prior(fm)
  expect_identical(dp$resp[dp$class == "lncor"], "Y")
  fp <- frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                 bf(w ~ x, family = gaussian())), data = dd,
            prior = set_prior("lkj(2)", class = "lncor", resp = "Y"))
  # an LKJ(2) pulls the correlation toward 0
  expect_lt(abs(summary(fp)$spec_pars["lncor_Y(b,c)", "Estimate"]),
            abs(rm[["Estimate"]]))
  expect_error(frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                        bf(w ~ x, family = gaussian())), data = dd,
                   prior = set_prior("lkj(2)", class = "lncor",
                                     resp = "w")),
               "No logistic-normal correlation matches")
  # predict() reads the response's own correlation block
  expect_identical(dim(predict(fm, resp = "Y")), c(300L, 4L, 3L))
})

test_that("two logistic_normal() responses keep their blocks apart", {
  dd <- smv_data(2)
  fm <- frm(mvbf(bf(Y ~ x, family = logistic_normal()),
                 bf(V ~ x, family = logistic_normal())), data = dd,
            control = smv_ctl)
  f1 <- frm(bf(Y ~ x), family = logistic_normal(), data = dd,
            control = smv_ctl)
  f2 <- frm(bf(V ~ x), family = logistic_normal(), data = dd,
            control = smv_ctl)
  ll <- as.numeric(logLik(f1)) + as.numeric(logLik(f2))
  expect_lt(abs(as.numeric(logLik(fm)) - ll), 1e-8 * abs(ll))
  sp <- summary(fm)$spec_pars
  expect_identical(grep("lncor", rownames(sp), value = TRUE),
                   c("lncor_Y(b,c)", "lncor_V(v,w)"))
  e2 <- summary(f2)$spec_pars["lncor(v,w)", ]
  expect_lt(abs(sp["lncor_V(v,w)", "Estimate"] - e2[["Estimate"]]),
            1e-3 * e2[["Est.Error"]])
  dp <- default_prior(fm)
  # in brms's row order, which sorts by response: brms 2.23.0 lists V's
  # lncor row before Y's (dev/brmsnames-1024-lncor.R)
  expect_identical(dp$resp[dp$class == "lncor"], c("V", "Y"))
})

test_that("dirichlet() in mvbf() is the sum of its two models", {
  dd <- smv_data(3)
  fm <- frm(mvbf(bf(V ~ x, family = dirichlet()),
                 bf(w ~ x, family = gaussian())), data = dd,
            control = smv_ctl)
  f1 <- frm(bf(V ~ x), family = dirichlet(), data = dd, control = smv_ctl)
  f2 <- frm(bf(w ~ x), data = dd, control = smv_ctl)
  ll <- as.numeric(logLik(f1)) + as.numeric(logLik(f2))
  expect_lt(abs(as.numeric(logLik(fm)) - ll), 1e-8 * abs(ll))
  # brms 2.23.0's fitted() of this model: 300 x 4 x 4, the composition's
  # P(Y = k) layers and then w (dev/catmulti-out/brms-fits.txt)
  fv <- fitted(fm)
  expect_identical(dim(fv), c(300L, 4L, 4L))
  expect_identical(dimnames(fv)[[3L]],
                   c("P(Y = u)", "P(Y = v)", "P(Y = w)", "w"))
  expect_lt(rc_rel_smv(fv[, "Estimate", 1:3], fitted(f1)[, "Estimate", ]),
            1e-5)
})


test_that("the correlation's name reads brms's Stan name of the response", {
  # Y_a is brms's Ya and cbind(a, b, c) its cbindabc; the first build
  # wrote lncor_Y_a and lncor_cbind(a, b, c), the second of which
  # hypothesis() could not parse (dev/reviews/2026-10-14-catmulti.md, 1)
  dd <- smv_data(4)
  dd$Y_a <- dd$Y
  f1 <- frm(mvbf(bf(Y_a ~ x, family = logistic_normal()),
                 bf(w ~ x, family = gaussian())), data = dd,
            control = smv_ctl)
  expect_true("lncor_Ya__b__c" %in% variables(f1))
  expect_true("lncor_Ya(b,c)" %in% rownames(summary(f1)$spec_pars))
  h1 <- hypothesis(f1, "lncor_Ya__b__c = 0", class = NULL)
  expect_true(is.finite(h1$hypothesis$Estimate))
  # set_prior() takes the response as frmtmb names it and as brms does
  for (r in c("Y_a", "Ya")) {
    fp <- frm(mvbf(bf(Y_a ~ x, family = logistic_normal()),
                   bf(w ~ x, family = gaussian())), data = dd,
              prior = set_prior("lkj(2)", class = "lncor", resp = r))
    expect_lt(abs(summary(fp)$spec_pars["lncor_Ya(b,c)", "Estimate"]),
              abs(summary(f1)$spec_pars["lncor_Ya(b,c)", "Estimate"]))
  }
  d2 <- cbind(dd, as.data.frame(unclass(dd$Y)))
  f2 <- frm(mvbf(bf(cbind(a, b, c) ~ x, family = logistic_normal()),
                 bf(w ~ x, family = gaussian())), data = d2,
            control = smv_ctl)
  expect_true("lncor_cbindabc__b__c" %in% variables(f2))
  h2 <- hypothesis(f2, "lncor_cbindabc__b__c = 0", class = NULL)
  expect_true(is.finite(h2$hypothesis$Estimate))
  expect_true(is.finite(h2$hypothesis$Est.Error))
})
