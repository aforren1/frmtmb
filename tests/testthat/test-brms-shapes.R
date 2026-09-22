# brms's return shapes on a maximum-likelihood fit (item 2.6f), and the
# predictive summary predict() became (item 2.6d).
#
# Every assertion here was run against frmtmb 0.60.0 first and recorded
# failing in dev/shapes-log/seen-failing-core.txt; the shapes were
# lme4's and glmmTMB's before this.

shapes_fit <- local({
  cache <- new.env(parent = emptyenv())
  function(key) {
    if (!is.null(cache[[key]])) return(cache[[key]])
    set.seed(20260917)
    n <- 120
    dd <- data.frame(x = rnorm(n), z = rnorm(n),
                     g = factor(rep(1:12, each = 10)))
    dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(12, 0, 0.7)[dd$g], 1)
    dd$cnt <- rpois(n, exp(0.4 + 0.3 * dd$x))
    dd$ord <- factor(cut(1 + 0.8 * dd$x + rnorm(n),
                         c(-Inf, -0.3, 0.8, Inf), labels = 1:3),
                     ordered = TRUE)
    out <- switch(key,
      gauss = frm(bf(y ~ x + (1 | g), sigma ~ z) + gaussian(), data = dd),
      pois = frm(bf(cnt ~ x) + poisson(), data = dd),
      ord = frm(bf(ord ~ x + cs(z)) + sratio(), data = dd),
      plain = frm(bf(y ~ x) + gaussian(), data = dd))
    cache[[key]] <- out
    out
  }
})

test_that("fixef() is brms's summary matrix, in brms's row order", {
  fit <- shapes_fit("gauss")
  fe <- fixef(fit)
  expect_true(is.matrix(fe))
  expect_equal(colnames(fe),
               c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  # brms puts every predictor's intercept first, then the ordinary
  # coefficients predictor by predictor
  expect_equal(rownames(fe),
               c("Intercept", "sigma_Intercept", "x", "sigma_z"))
  expect_equal(unname(fe[, "Estimate"]),
               unname(fixef(fit, flatten = TRUE)[
                 c("(Intercept)", "sigma_(Intercept)", "x", "sigma_z")]))
  # Est.Error is the standard error and the Q columns a Wald interval
  expect_equal(unname(fe[, "Est.Error"]),
               unname(sqrt(diag(vcov(fit)))[rownames(fe)]))
  expect_equal(unname(fe[, "Q2.5"]),
               unname(fe[, "Estimate"] - stats::qnorm(0.975) *
                        fe[, "Est.Error"]))
})

test_that("fixef() takes brms's pars and probs", {
  fit <- shapes_fit("gauss")
  expect_equal(rownames(fixef(fit, pars = "x")), "x")
  expect_equal(colnames(fixef(fit, probs = c(0.1, 0.9))),
               c("Estimate", "Est.Error", "Q10", "Q90"))
  expect_error(fixef(fit, pars = "nope"), "no such parameter")
})

test_that("fixef() still refuses the draws-only arguments by name", {
  fit <- shapes_fit("gauss")
  expect_error(fixef(fit, summary = FALSE), "summary = FALSE")
  expect_error(fixef(fit, robust = TRUE), "robust = TRUE")
})

test_that("vcov() covers brms's population-level coefficients only", {
  # an intercept-only dpar nobody wrote a formula for is brms's
  # spec_par, not a population-level coefficient: it is why vcov() was
  # one row wider than brms's on the port's fixture
  fit <- shapes_fit("plain")
  expect_equal(rownames(vcov(fit)), c("Intercept", "x"))
  expect_false("sigma" %in% rownames(vcov(fit)))
  fit2 <- shapes_fit("gauss")
  expect_equal(rownames(vcov(fit2)),
               c("Intercept", "sigma_Intercept", "x", "sigma_z"))
})

test_that("vcov() takes brms's correlation and pars", {
  fit <- shapes_fit("plain")
  cr <- vcov(fit, correlation = TRUE)
  expect_equal(unname(diag(cr)), rep(1, nrow(cr)))
  # brms spells it `correlation`, so the partial match `cor` reaches it
  expect_equal(vcov(fit, cor = TRUE), cr)
  expect_equal(dim(vcov(fit, pars = "x")), c(1L, 1L))
})

test_that("ngrps() is brms's named list", {
  fit <- shapes_fit("gauss")
  expect_equal(ngrps(fit), list(g = 12L))
  # brms has no group-level entry for a fit without one
  expect_null(ngrps(shapes_fit("plain")))
})

test_that("fitted() is brms's four-column summary", {
  fit <- shapes_fit("gauss")
  fi <- fitted(fit)
  expect_equal(dim(fi), c(nobs(fit), 4L))
  expect_equal(colnames(fi), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  # the Estimate is the value fitted() reported before the shape moved
  expect_equal(unname(fi[, "Estimate"]),
               unname(frm_linpred(fit, type = "response")))
  expect_equal(unname(fi[, "Est.Error"]),
               unname(frm_linpred(fit, type = "response",
                                  se.fit = TRUE)$se.fit))
  expect_equal(colnames(fitted(fit, probs = c(0.1, 0.9))),
               c("Estimate", "Est.Error", "Q10", "Q90"))
})

test_that("fitted() on an ordinal fit is brms's n x 4 x K array", {
  fit <- shapes_fit("ord")
  fi <- fitted(fit)
  expect_equal(dim(fi), c(nobs(fit), 4L, 3L))
  expect_equal(dimnames(fi)[[2]],
               c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_equal(dimnames(fi)[[3]], paste0("P(Y = ", 1:3, ")"))
  # the estimates still sum to one across the categories
  expect_equal(unname(rowSums(fi[, "Estimate", ])), rep(1, nobs(fit)))
  # and the category probabilities carry a standard error now
  expect_true(all(is.finite(fi[, "Est.Error", ])))
  expect_true(all(fi[, "Est.Error", ] > 0))
})

test_that("fitted() refuses the arguments that need draws", {
  fit <- shapes_fit("gauss")
  expect_error(fitted(fit, ndraws = 3), "ndraws")
  expect_error(fitted(fit, summary = FALSE), "summary")
  expect_error(fitted(fit, robust = TRUE), "robust")
})

test_that("fitted() takes nlpar and allow_new_levels", {
  fit <- shapes_fit("gauss")
  nd <- data.frame(x = 0, z = 0, g = factor("new"))
  expect_error(fitted(fit, newdata = nd), "New levels")
  fi <- fitted(fit, newdata = nd, allow_new_levels = TRUE)
  expect_equal(dim(fi), c(1L, 4L))
  # nlpar is brms's spelling of a non-linear parameter, which this
  # package reaches through dpar
  expect_equal(fitted(fit, nlpar = "sigma"), fitted(fit, dpar = "sigma"))
})

test_that("residuals() is brms's four-column summary", {
  fit <- shapes_fit("gauss")
  rs <- residuals(fit)
  expect_equal(dim(rs), c(nobs(fit), 4L))
  expect_equal(colnames(rs), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  # the estimate is the raw residual: the response minus the fitted
  # value, which is what residuals() returned before the shape moved
  expect_equal(unname(rs[, "Estimate"]),
               unname(model.frame(fit)[[1L]] - fitted(fit)[, "Estimate"]))
  expect_equal(dim(residuals(fit, type = "pearson", probs = 0.65)),
               c(nobs(fit), 3L))
})

test_that("summary() carries brms's fixed and random slots", {
  fit <- shapes_fit("gauss")
  s <- summary(fit)
  expect_true(is.data.frame(s$fixed))
  expect_equal(rownames(s$fixed),
               c("Intercept", "sigma_Intercept", "x", "sigma_z"))
  # brms's four columns first, then the Wald test brms has no
  # counterpart for (it writes Rhat and the two ESS columns there)
  expect_equal(colnames(s$fixed),
               c("Estimate", "Est.Error", "l-95% CI", "u-95% CI",
                 "z value", "Pr(>|z|)"))
  expect_true(is.list(s$random))
  expect_equal(rownames(s$random$g), "sd(Intercept)")
  expect_equal(s$nobs, nobs(fit))
  expect_equal(s$ngrps, list(g = 12L))
  expect_output(print(s), "Regression Coefficients:")
  expect_output(print(s), "Multilevel Hyperparameters:")
})

test_that("summary(priors = TRUE) prints the priors, as brms does", {
  fit <- shapes_fit("plain")
  s <- summary(fit, priors = TRUE)
  expect_output(print(s), "Priors:")
  expect_false(any(grepl("Priors:", capture.output(print(summary(fit))))))
})

test_that("print() of a fit uses brms's section headings", {
  fit <- shapes_fit("gauss")
  expect_output(print(fit), "Multilevel Hyperparameters:")
  expect_output(print(fit), "Regression Coefficients:")
})

test_that("variables() lists an ordinal fit's thresholds and cs terms", {
  fit <- shapes_fit("ord")
  v <- variables(fit)
  expect_true(all(c("b_Intercept[1]", "b_Intercept[2]") %in% v))
  expect_true(all(c("bcs_z[1]", "bcs_z[2]") %in% v))
  # the thresholds reported are the model's own and not the internal
  # (first threshold, log increment) parameterization sratio estimates,
  # so they are increasing
  th <- vapply(c("b_Intercept[1]", "b_Intercept[2]"),
               function(nm) hypothesis(fit, paste0("`", nm, "` = 0"),
                                       class = NULL)$hypothesis$Estimate,
               numeric(1))
  expect_true(th[2] > th[1])
  expect_false(isTRUE(all.equal(unname(th),
                                unname(fit$estimates[["tau_raw"]][1:2]))))
})

test_that("predict() is brms's predictive summary", {
  fit <- shapes_fit("pois")
  set.seed(1)
  pr <- predict(fit)
  expect_equal(dim(pr), c(nobs(fit), 4L))
  expect_equal(colnames(pr), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  # observation noise is in it: the predictive spread is far wider than
  # the standard error of the mean
  fi <- fitted(fit)
  expect_true(mean(pr[, "Est.Error"]) > 3 * mean(fi[, "Est.Error"]))
  # and a count predictive interval is on the counts
  expect_true(all(pr[, "Q2.5"] >= 0))
  expect_equal(dim(predict(fit, probs = c(0.2, 0.5, 0.8))),
               c(nobs(fit), 5L))
})

test_that("predict(summary = FALSE) gives the simulated draws", {
  fit <- shapes_fit("pois")
  set.seed(2)
  d <- predict(fit, summary = FALSE, ndraws = 25)
  expect_equal(dim(d), c(25L, nobs(fit)))
})

test_that("predict() on an ordinal fit gives brms's P(Y = k) columns", {
  fit <- shapes_fit("ord")
  set.seed(3)
  pr <- predict(fit)
  expect_equal(dim(pr), c(nobs(fit), 3L))
  expect_equal(colnames(pr), paste0("P(Y = ", 1:3, ")"))
  expect_equal(unname(rowSums(pr)), rep(1, nobs(fit)))
})

test_that("predict() refuses the retired glmmTMB type vocabulary", {
  fit <- shapes_fit("pois")
  expect_error(predict(fit, type = "link"), "frm_linpred")
  expect_error(predict(fit, se.fit = TRUE), "frm_linpred")
})

test_that("frm_linpred() is the old predict(), unchanged", {
  fit <- shapes_fit("gauss")
  expect_equal(frm_linpred(fit, type = "response"),
               unname(fitted(fit)[, "Estimate"]),
               ignore_attr = TRUE)
  p <- frm_linpred(fit, se.fit = TRUE)
  expect_named(p, c("fit", "se.fit"))
  expect_equal(frm_linpred(fit, type = "link"),
               frm_linpred(fit, dpar = "mu"))
})

test_that("predict() covers at close to the nominal rate", {
  # the whole point of simulating from the joint precision rather than
  # from the plug-in estimates alone; the number of replicates here is
  # a smoke test, and dev/shapes-coverage.R is the measurement
  set.seed(20260917)
  hit <- 0L
  tot <- 0L
  for (r in 1:20) {
    n <- 60
    dd <- data.frame(x = rnorm(n))
    dd$y <- rnorm(n, 1 + 0.5 * dd$x, 1)
    f <- frm(bf(y ~ x) + gaussian(), data = dd)
    nd <- data.frame(x = rnorm(10))
    pr <- predict(f, newdata = nd, ndraws = 400)
    ynew <- rnorm(10, 1 + 0.5 * nd$x, 1)
    hit <- hit + sum(ynew >= pr[, "Q2.5"] & ynew <= pr[, "Q97.5"])
    tot <- tot + 10L
  }
  # 200 predictions at 0.95: the 0.999 binomial lower bound is 0.90
  expect_gt(hit / tot, 0.90)
})
