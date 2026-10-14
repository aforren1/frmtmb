# brms 2.23.0's simplex-response families: dirichlet(), logistic_normal()
# and brmsfamily("dirichlet2").
#
# The reference densities here are written from lgamma(), log() and a
# multivariate normal density by hand, after brms's Stan functions
# (chunks/fun_dirichlet_logit.stan, fun_logistic_normal.stan, and Stan's
# own dirichlet_lpdf for dirichlet2), not from the package code. The
# brms behaviour each test pins was read off brms 2.23.0 itself:
# dev/famsimplex-brms-probe.R and dev/famsimplex-brms-fits2.R, recorded
# in dev/famsimplex-findings.md. Tolerances are ratios: to the double
# epsilon times the sum of the absolute terms for a density (the
# lgamma() terms cancel, so the value alone does not measure what
# rounding can reach), to the measured log-likelihood for a fit, and to
# the measured standard error for an estimate.

ULP <- .Machine$double.eps

# the error of a density against its reference, in units of the double
# epsilon times the reference's own sum of absolute terms
dens_ulps <- function(ours, ref) {
  max(abs(ours - as.numeric(ref)) / attr(ref, "scale")) / ULP
}

sx_softmax <- function(E) {
  E <- E - apply(E, 1L, max)
  P <- exp(E)
  P / rowSums(P)
}

sx_rdirichlet <- function(A) {
  G <- matrix(stats::rgamma(length(A), shape = A), nrow(A))
  G / rowSums(G)
}

ref_dirichlet <- function(Y, A) {
  t1 <- lgamma(rowSums(A))
  t2 <- lgamma(A)
  t3 <- (A - 1) * log(Y)
  structure(t1 - rowSums(t2) + rowSums(t3),
            scale = abs(t1) + rowSums(abs(t2)) + rowSums(abs(t3)))
}

ref_logistic_normal <- function(Y, M, S, C, ref) {
  X <- log(Y[, -ref, drop = FALSE]) - log(Y[, ref])
  D <- ncol(X)
  parts <- vapply(seq_len(nrow(Y)), function(i) {
    Sig <- diag(S[i, ], D) %*% C %*% diag(S[i, ], D)
    r <- X[i, ] - M[i, ]
    c(D * log(2 * pi), as.numeric(determinant(Sig)$modulus),
      sum(r * solve(Sig, r)), sum(log(Y[i, ])), sum(abs(log(Y[i, ]))))
  }, numeric(5))
  structure(-0.5 * colSums(parts[1:3, , drop = FALSE]) - parts[4, ],
            scale = 0.5 * colSums(abs(parts[1:3, , drop = FALSE])) +
              parts[5, ])
}

sx_data_dir <- function(seed, n = 300, ng = 0) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  E <- cbind(0, 0.4 + 0.6 * dd$x, -0.3 - 0.5 * dd$x)
  if (ng > 0) {
    dd$g <- factor(rep(seq_len(ng), length.out = n))
    U <- matrix(stats::rnorm(ng * 2), ng) %*%
      chol(matrix(c(0.25, 0.06, 0.06, 0.09), 2))
    E[, 2:3] <- E[, 2:3] + U[as.integer(dd$g), ]
  }
  Y <- sx_rdirichlet(sx_softmax(E) * exp(2 + 0.5 * dd$z))
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd
}

sx_data_ln <- function(seed, n = 300, lncor = 0.3) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  M <- cbind(0.4 + 0.6 * dd$x, -0.3 - 0.5 * dd$x)
  sb <- exp(-0.5 + 0.4 * dd$z)
  sc <- exp(-0.2)
  W <- matrix(stats::rnorm(n * 2), n) %*%
    chol(matrix(c(1, lncor, lncor, 1), 2))
  Y <- sx_softmax(cbind(0, M[, 1] + sb * W[, 1], M[, 2] + sc * W[, 2]))
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd
}

# ---------------------------------------------------------------------
# The family objects
# ---------------------------------------------------------------------

test_that("the constructors answer brms's link fields", {
  # brms 2.23.0: dirichlet()$link "logit", $link_phi "log";
  # logistic_normal()$link "identity", $link_sigma "log";
  # brmsfamily("dirichlet2")$link "log"
  expect_identical(dirichlet()$link, "logit")
  expect_identical(dirichlet()$link_phi, "log")
  expect_identical(dirichlet(link_phi = "softplus")$link_phi, "softplus")
  expect_identical(logistic_normal()$link, "identity")
  expect_identical(logistic_normal()$link_sigma, "log")
  expect_identical(brmsfamily("dirichlet2")$link, "log")
  expect_identical(brmsfamily("dirichlet2", "softplus")$link, "softplus")
  expect_identical(brmsfamily("dirichlet")$family, "dirichlet")
  expect_identical(brmsfamily("logistic_normal")$family, "logistic_normal")
})

test_that("links and refcat are refused in brms's words", {
  expect_error(dirichlet(link = "probit"),
               "'probit' is not a supported link for family 'dirichlet'")
  expect_error(logistic_normal(link = "logit"),
               "not a supported link for family 'logistic_normal'")
  expect_error(dirichlet(link_phi = "logit"),
               "'logit' is not a supported link for parameter 'phi'")
  expect_error(brmsfamily("dirichlet2", "logit"),
               "not a supported link for family 'dirichlet2'")
  # brms allows refcat = NA for dirichlet() and refuses it for
  # logistic_normal(), whose log ratios need a reference
  expect_error(logistic_normal(refcat = NA), "needs a reference category")
  expect_error(dirichlet(refcat = c("a", "b")), "with one string")
  expect_s3_class(dirichlet(refcat = NA), "frmtmb_family")
  dd <- sx_data_dir(1)
  expect_error(frm(bf(Y ~ x), family = dirichlet(refcat = "zz"), data = dd),
               "The reference response category must be one of a, b, c")
  # brmsfamily() takes refcat for these two families now, and still
  # refuses it on a family that has none
  expect_error(brmsfamily("gaussian", refcat = "a"), "has no argument")
})

# ---------------------------------------------------------------------
# The response
# ---------------------------------------------------------------------

test_that("the response check follows brms's, row by row", {
  dd <- sx_data_dir(2)
  Y <- dd$Y
  fit_y <- function(Yb, fam = dirichlet()) {
    d <- dd
    d$Y <- Yb
    frm(bf(Y ~ x), family = fam, data = d)
  }
  # the ABSENT case: a composition, and one rounded to 1e-9 of 1, fit
  expect_s3_class(fit_y(Y), "frmtmb_fit")
  Yb <- Y
  Yb[1, 1] <- Yb[1, 1] + 1e-9
  expect_s3_class(fit_y(Yb), "frmtmb_fit")
  # one row 1e-5 off the simplex: brms refuses it too ("must sum to 1")
  Yb <- Y
  Yb[1, 1] <- Yb[1, 1] + 1e-5
  expect_error(fit_y(Yb), "Response values in simplex models must sum to 1")
  # ONE row 1e-7 off in 300: brms's all.equal() test is a mean over rows
  # and passes it (dev/famsimplex-edge.log); the per-row test refuses it
  Yb <- Y
  Yb[1, 1] <- Yb[1, 1] + 1e-7
  expect_error(fit_y(Yb, logistic_normal()),
               "1 row does not [(]row 1;")
  # zeros and negatives: brms's ybounds message
  Yb <- Y
  Yb[1, ] <- c(0.5, 0.5, 0)
  expect_error(fit_y(Yb), "Family 'dirichlet' requires response greater than 0")
  Yb[1, ] <- c(1.2, -0.2, 0)
  expect_error(fit_y(Yb, logistic_normal()),
               "Family 'logistic_normal' requires response greater than 0")
  # a row with an NA is dropped, as brms drops it
  Yb <- Y
  Yb[1, 1] <- NA
  f <- suppressMessages(fit_y(Yb))
  expect_identical(nobs(f), nrow(Y) - 1L)
  # one category, or a vector: brms's "At least 2 response categories"
  expect_error(fit_y(Y[, 1, drop = FALSE]),
               "At least 2 response categories are required")
  d2 <- dd
  d2$y1 <- Y[, 1]
  expect_error(frm(bf(y1 ~ x), family = dirichlet(), data = d2),
               "At least 2 response categories are required")
  # only weights() is read, as in brms
  d2$w <- 2
  expect_error(frm(bf(Y | trials(w) ~ x), family = dirichlet(), data = d2),
               "trials\\(\\)")
})

test_that("the categories come from the response's columns, as in brms", {
  dd <- sx_data_dir(3)
  # cbind(y1, y2, y3) names them y1, y2, y3, which brms turns into the
  # dpar names muy2 and muy3; an unnamed matrix gives mu2 and mu3
  dd$y1 <- dd$Y[, 1]
  dd$y2 <- dd$Y[, 2]
  dd$y3 <- dd$Y[, 3]
  f1 <- frm(bf(cbind(y1, y2, y3) ~ x), family = dirichlet(), data = dd)
  expect_identical(rownames(fixef(f1)),
                   c("muy2_Intercept", "muy3_Intercept", "muy2_x",
                     "muy3_x"))
  d2 <- dd
  d2$Y <- unname(dd$Y)
  f2 <- frm(bf(Y ~ x), family = dirichlet(), data = d2)
  expect_identical(rownames(fixef(f2))[1:2],
                   c("mu2_Intercept", "mu3_Intercept"))
  # the same model under three spellings of its response
  expect_lt(abs(as.numeric(logLik(f1) - logLik(f2))),
            1e-8 * abs(as.numeric(logLik(f2))))
})

# ---------------------------------------------------------------------
# The densities, against references written by hand
# ---------------------------------------------------------------------

test_that("the dirichlet density is brms's dirichlet_logit_lpdf()", {
  set.seed(4)
  for (K in 2:4) {
    cats <- letters[seq_len(K)]
    for (ref in c(1L, K, NA_integer_)) {
      fam <- frmtmb:::fam_dirichlet_impl(cats, ref)
      sx <- fam$simplex
      n <- 50
      Y <- sx_rdirichlet(matrix(1.5, n, K))
      E <- matrix(stats::rnorm(n * length(sx$pred), 0, 1.5), n)
      phi <- exp(stats::runif(n, -1, 4))
      M <- matrix(0, n, K)
      M[, sx$pred] <- E
      dp <- c(stats::setNames(lapply(seq_along(sx$pred),
                                     function(j) E[, j]), sx$mu),
              list(phi = phi))
      ours <- fam$lpdf(Y, dp, list())
      ref_v <- ref_dirichlet(Y, sx_softmax(M) * phi)
      expect_lt(dens_ulps(ours, ref_v), 64,
                label = paste0("dirichlet K = ", K, " ref = ", ref))
    }
  }
})

test_that("the dirichlet2 density is Stan's dirichlet_lpdf() on alpha", {
  set.seed(5)
  fam <- frmtmb:::fam_dirichlet2_impl(c("a", "b", "c"))
  Y <- sx_rdirichlet(matrix(2, 60, 3))
  A <- matrix(exp(stats::runif(180, -2, 4)), 60)
  dp <- stats::setNames(lapply(1:3, function(k) A[, k]), fam$simplex$mu)
  ours <- fam$lpdf(Y, dp, list())
  ref_v <- ref_dirichlet(Y, A)
  expect_lt(dens_ulps(ours, ref_v), 64)
})

test_that("the logistic-normal density is brms's, Jacobian included", {
  set.seed(6)
  for (K in 2:4) {
    cats <- letters[seq_len(K)]
    D <- K - 1L
    for (ref in unique(c(1L, K))) {
      fam <- frmtmb:::fam_logistic_normal_impl(cats, ref)
      sx <- fam$simplex
      n <- 40
      Y <- sx_rdirichlet(matrix(2, n, K))
      M <- matrix(stats::rnorm(n * D), n)
      S <- matrix(exp(stats::runif(n * D, -1, 1)), n)
      C <- if (D > 1L) {
        stats::cov2cor(crossprod(matrix(stats::rnorm(D * (D + 3)), D + 3)))
      } else matrix(1, 1, 1)
      dp <- c(stats::setNames(lapply(seq_len(D), function(j) M[, j]), sx$mu),
              stats::setNames(lapply(seq_len(D), function(j) S[, j]),
                              sx$sigma))
      ex <- if (D > 1L) list(lncor_raw = frmtmb:::us_theta_cor(C)) else NULL
      ours <- if (is.null(ex)) fam$lpdf(Y, dp, list()) else {
        fam$lpdf(Y, dp, list(), ex)
      }
      ref_v <- ref_logistic_normal(Y, M, S, C, ref)
      expect_lt(dens_ulps(ours, ref_v), 64,
                label = paste0("logistic_normal K = ", K, " ref = ", ref))
    }
  }
})

# ---------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------

test_that("a dirichlet fit's logLik is the hand density at its estimates", {
  dd <- sx_data_dir(7)
  fit <- frm(bf(Y ~ x, phi ~ z), family = dirichlet(), data = dd)
  b <- fixef(fit)[, "Estimate"]
  E <- cbind(0, b[["mub_Intercept"]] + b[["mub_x"]] * dd$x,
             b[["muc_Intercept"]] + b[["muc_x"]] * dd$x)
  phi <- exp(b[["phi_Intercept"]] + b[["phi_z"]] * dd$z)
  ll <- sum(ref_dirichlet(dd$Y, sx_softmax(E) * phi))
  expect_lt(abs(ll - as.numeric(logLik(fit))), 1e-8 * abs(ll))
  # the truth is inside four standard errors of every estimate
  truth <- c(mub_Intercept = 0.4, muc_Intercept = -0.3, phi_Intercept = 2,
             mub_x = 0.6, muc_x = -0.5, phi_z = 0.5)
  z <- (b[names(truth)] - truth) / fixef(fit)[names(truth), "Est.Error"]
  expect_lt(max(abs(z)), 4)
  # brms's Links line: mub = logit, as for categorical in brms
  expect_match(paste(utils::capture.output(print(fit)), collapse = "\n"),
               "mub = logit; muc = logit; phi = log", fixed = TRUE)
})

test_that("changing the reference category reparameterizes, as in brms", {
  dd <- sx_data_dir(8)
  f1 <- frm(bf(Y ~ x), family = dirichlet(), data = dd)
  f2 <- frm(bf(Y ~ x), family = dirichlet(refcat = "b"), data = dd)
  expect_identical(rownames(fixef(f2))[1:2],
                   c("mua_Intercept", "muc_Intercept"))
  ll <- abs(as.numeric(logLik(f1)))
  expect_lt(abs(as.numeric(logLik(f1) - logLik(f2))), 1e-8 * ll)
  # mua under reference b is minus mub under reference a
  e1 <- fixef(f1)[, "Estimate"]
  e2 <- fixef(f2)[, "Estimate"]
  expect_lt(abs(e2[["mua_x"]] + e1[["mub_x"]]),
            1e-3 * fixef(f1)["mub_x", "Est.Error"])
  # a brms family object carries its refcat through
  skip_if_not_installed("brms")
  f3 <- frm(bf(Y ~ x), family = brms::dirichlet(refcat = "b"), data = dd)
  expect_identical(rownames(fixef(f3)), rownames(fixef(f2)))
})

test_that("refcat = NA leaves a common shift unidentified, and says so", {
  # brms allows it and leans on its priors; at the maximum the
  # likelihood is flat along the shift, which the SE check names
  dd <- sx_data_dir(9)
  w <- character(0)
  f <- withCallingHandlers(
    frm(bf(Y ~ x), family = dirichlet(refcat = NA), data = dd),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("Standard errors are not available", w)))
  expect_identical(rownames(fixef(f))[1:3],
                   c("mua_Intercept", "mub_Intercept", "muc_Intercept"))
  # the fitted means are still the reference-a fit's, to the optimizer
  f1 <- frm(bf(Y ~ x), family = dirichlet(), data = dd)
  expect_lt(abs(as.numeric(logLik(f) - logLik(f1))),
            1e-6 * abs(as.numeric(logLik(f1))))
})

test_that("logistic_normal: lncor is reported as brms reports it", {
  dd <- sx_data_ln(10)
  fit <- frm(bf(Y ~ x, sigmab ~ z), family = logistic_normal(), data = dd)
  sp <- summary(fit)$spec_pars
  expect_identical(rownames(sp), c("sigmac", "lncor(b,c)"))
  # the delta-method error and the Fisher-z interval bracket the estimate
  r <- sp["lncor(b,c)", ]
  expect_true(r[["l-95% CI"]] < r[["Estimate"]] &&
                r[["Estimate"]] < r[["u-95% CI"]])
  expect_lt(abs(r[["Estimate"]] - 0.3), 4 * r[["Est.Error"]])
  # hypothesis() reads it under brms's variable name
  h <- hypothesis(fit, "lncor__b__c > 0", class = NULL)
  expect_lt(abs(h$hypothesis$Estimate - r[["Estimate"]]),
            1e-8 * abs(r[["Estimate"]]))
  # the maximum does not depend on the reference category: the
  # covariance is unrestricted, so another reference is a linear change
  # of variables with the same Jacobian on the simplex
  f2 <- frm(bf(Y ~ x, sigmab ~ z), family = logistic_normal(refcat = "c"),
            data = dd)
  expect_identical(rownames(summary(f2)$spec_pars)[2], "lncor(a,b)")
  f1 <- frm(bf(Y ~ x), family = logistic_normal(), data = dd)
  f3 <- frm(bf(Y ~ x), family = logistic_normal(refcat = "c"), data = dd)
  expect_lt(abs(as.numeric(logLik(f1) - logLik(f3))),
            1e-6 * abs(as.numeric(logLik(f1))))
})

test_that("group effects correlate across categories with (1 | ID | g)", {
  dd <- sx_data_dir(11, n = 400, ng = 20)
  fit <- frm(bf(Y ~ x + (1 | ID | g)), family = dirichlet(), data = dd)
  vc <- VarCorr(fit)[["g"]]
  expect_identical(rownames(vc$sd), c("mub_Intercept", "muc_Intercept"))
  expect_true(is.finite(vc$cor["muc_Intercept", "Estimate",
                               "mub_Intercept"]))
  # without the ID the two categories' effects are separate blocks
  f0 <- frm(bf(Y ~ x + (1 | g)), family = dirichlet(), data = dd)
  expect_null(dim(VarCorr(f0)[["g"]]$cor))
  f <- fitted(fit)
  expect_identical(dim(f), c(400L, 4L, 3L))
})

# ---------------------------------------------------------------------
# Post-fit routes, each in brms's shape
# ---------------------------------------------------------------------

test_that("fitted() and predict() return brms's arrays", {
  dd <- sx_data_dir(12)
  fit <- frm(bf(Y ~ x, phi ~ z), family = dirichlet(), data = dd)
  f <- fitted(fit)
  # brms 2.23.0: 300 x 4 x 3, P(Y = a) ... (dev/famsimplex-brms-fits2-dir.log)
  expect_identical(dim(f), c(300L, 4L, 3L))
  expect_identical(dimnames(f)[[3]], c("P(Y = a)", "P(Y = b)", "P(Y = c)"))
  expect_null(dimnames(f)[[1]])
  expect_lt(max(abs(rowSums(f[, "Estimate", ]) - 1)), 64 * ULP)
  expect_true(all(f[, "Est.Error", ] > 0))
  p <- predict(fit)
  expect_identical(dim(p), c(300L, 4L, 3L))
  expect_identical(dimnames(p)[[3]], c("a", "b", "c"))
  # simulated compositions: every row on the simplex
  s <- simulate(fit, nsim = 2, seed = 1)
  Ys <- s[[1]]
  expect_identical(colnames(Ys), c("a", "b", "c"))
  expect_lt(max(abs(rowSums(Ys) - 1)), 64 * ULP)
  expect_true(all(Ys > 0))
  # the latent predictor of one category
  lp <- frm_linpred(fit, dpar = "mub")
  b <- fixef(fit)[, "Estimate"]
  expect_lt(max(abs(lp - (b[["mub_Intercept"]] + b[["mub_x"]] * dd$x))),
            1e-10 * max(abs(lp)))
})

test_that("fitted() refuses logistic_normal in brms's words", {
  dd <- sx_data_ln(13)
  fit <- frm(bf(Y ~ x), family = logistic_normal(), data = dd)
  msg <- paste0("Cannot compute expected values of the posterior ",
                "predictive distribution for family 'logistic_normal'")
  expect_error(fitted(fit), msg, fixed = TRUE)
  expect_error(conditional_effects(fit, "x"), msg, fixed = TRUE)
  # what does work: draws on the simplex and the category predictors
  p <- predict(fit)
  expect_identical(dimnames(p)[[3]], c("a", "b", "c"))
  Ys <- simulate(fit, nsim = 1, seed = 2)[[1]]
  expect_lt(max(abs(rowSums(Ys) - 1)), 64 * ULP)
})

test_that("residuals() and pp_check() are refused, as brms refuses them", {
  dd <- sx_data_dir(14)
  fit <- frm(bf(Y ~ x), family = dirichlet(), data = dd)
  for (ty in c("response", "pearson", "deviance")) {
    expect_error(residuals(fit, type = ty), "not defined for")
  }
  skip_if_not_installed("bayesplot")
  expect_error(pp_check(fit), "'pp_check' is not implemented for this family")
})

test_that("conditional_effects() draws one curve per category", {
  dd <- sx_data_dir(15)
  fit <- frm(bf(Y ~ x), family = dirichlet(), data = dd)
  ce <- conditional_effects(fit, "x")[[1]]
  expect_identical(levels(ce$cats__), c("a", "b", "c"))
  # the curve is fitted() on the grid, Estimate and Est.Error alike
  nd <- ce[ce$cats__ == "b", , drop = FALSE]
  f <- fitted(fit, newdata = nd, re_formula = NA)
  expect_lt(max(abs(nd$estimate__ - f[, "Estimate", 2])), 64 * ULP)
  expect_lt(max(abs(nd$se__ - f[, "Est.Error", 2])), 64 * ULP)
  expect_true(all(nd$lower__ < nd$estimate__ & nd$estimate__ < nd$upper__))
  expect_error(conditional_effects(fit, "x", categorical = FALSE),
               "nothing to draw")
})

test_that("the prior table has brms's classes, and lncor takes lkj()", {
  dd <- sx_data_ln(16)
  tab <- default_prior(bf(Y ~ x, sigmab ~ z), data = dd,
                       family = logistic_normal())
  expect_true("lncor" %in% tab$class)
  expect_true("sigmac" %in% tab$class)
  expect_true(all(c("mub", "muc", "sigmab") %in% tab$dpar))
  # a MAP fit under lkj(4) pulls the correlation toward zero
  f0 <- frm(bf(Y ~ x), family = logistic_normal(), data = dd)
  f4 <- frm(bf(Y ~ x), family = logistic_normal(), data = dd,
            prior = set_prior("lkj(4)", class = "lncor"))
  r0 <- summary(f0)$spec_pars["lncor(b,c)", "Estimate"]
  r4 <- summary(f4)$spec_pars["lncor(b,c)", "Estimate"]
  expect_lt(abs(r4), abs(r0))
  expect_error(set_prior("normal(0, 1)", class = "lncor"),
               "takes an lkj[(][)] prior")
  # a model without the correlations refuses the class
  dd2 <- sx_data_dir(16)
  expect_error(frm(bf(Y ~ x), family = dirichlet(), data = dd2,
                   prior = set_prior("lkj(2)", class = "lncor")),
               "No logistic-normal correlation")
})

test_that("the first predicted category may take its own formula", {
  # brms 2.23.0 builds bf(Y ~ x, mub ~ z) on a dirichlet() response
  # whose first predicted category is b (dev/famsimplex-brms-first.log);
  # the parser used to keep the first location dpar for the main formula
  # alone, on every family with several of them
  dd <- sx_data_dir(17)
  fit <- frm(bf(Y ~ x, mub ~ z), family = dirichlet(), data = dd)
  expect_identical(rownames(fixef(fit)),
                   c("mub_Intercept", "muc_Intercept", "mub_z", "muc_x"))
})

test_that("categorical() takes a formula for its first predicted category", {
  # the same model under a level order that makes b the SECOND predicted
  # category, which was always allowed; brms 2.23.0 itself stops on this
  # call inside use_glm_primitive_categorical() (findings, brms bugs)
  set.seed(18)
  n <- 300
  d2 <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  E <- cbind(0, 0.3 + 0.5 * d2$z, -0.2 + 0.6 * d2$x)
  P <- sx_softmax(E)
  d2$f <- factor(c("a", "b", "c")[apply(P, 1L, function(p) {
    sample.int(3L, 1L, prob = p)
  })], levels = c("a", "b", "c"))
  f1 <- frm(bf(f ~ x, mub ~ z), family = categorical(), data = d2)
  d3 <- d2
  d3$f <- factor(d3$f, levels = c("a", "c", "b"))
  f2 <- frm(bf(f ~ x, mub ~ z), family = categorical(), data = d3)
  expect_lt(abs(as.numeric(logLik(f1) - logLik(f2))),
            1e-6 * abs(as.numeric(logLik(f2))))
  expect_lt(abs(fixef(f1)["mub_z", "Estimate"] -
                  fixef(f2)["mub_z", "Estimate"]),
            1e-3 * fixef(f2)["mub_z", "Est.Error"])
})

test_that("a constant on a location of a family with several is refused", {
  # the location was built from a formula whatever the constant said, so
  # the fit returned the estimate the user had asked to replace; brms
  # 2.23.0 ignores it the same way (dev/famsimplex-rev-pfix-*.log)
  held <- "cannot be held at a value"
  dd <- sx_data_dir(19)
  expect_error(frm(bf(Y ~ x, muc = 0), family = dirichlet(), data = dd),
               held)
  expect_error(frm(bf(Y ~ x, mub = 0), family = dirichlet(), data = dd),
               "mub = 0")
  dl <- sx_data_ln(19)
  expect_error(frm(bf(Y ~ x, muc = 0), family = logistic_normal(),
                   data = dl), held)
  # the ABSENT case: a constant on a dpar that is not a location is
  # held, and the fit reads it
  f5 <- frm(bf(Y ~ x, phi = 5), family = dirichlet(), data = dd)
  ff <- frm(bf(Y ~ x), family = dirichlet(), data = dd)
  expect_lt(as.numeric(logLik(f5)), as.numeric(logLik(ff)))
  E <- cbind(0, fixef(f5)["mub_Intercept", "Estimate"] +
               fixef(f5)["mub_x", "Estimate"] * dd$x,
             fixef(f5)["muc_Intercept", "Estimate"] +
               fixef(f5)["muc_x", "Estimate"] * dd$x)
  ll5 <- sum(ref_dirichlet(dd$Y, sx_softmax(E) * 5))
  expect_lt(abs(ll5 - as.numeric(logLik(f5))), 1e-8 * abs(ll5))
  expect_s3_class(frm(bf(Y ~ x, sigmac = 1), family = logistic_normal(),
                      data = dl), "frmtmb_fit")
  # the refusal names the route that does hold a location: an
  # offset-only formula, which the fit then reads exactly
  expect_error(frm(bf(Y ~ x, muc = 0.2), family = dirichlet(), data = dd),
               "muc ~ 0 + offset(o)", fixed = TRUE)
  dd$o <- 0.2
  fo <- frm(bf(Y ~ x, muc ~ 0 + offset(o)), family = dirichlet(), data = dd)
  expect_false(any(startsWith(rownames(fixef(fo)), "muc_")))
  b <- fixef(fo)[, "Estimate"]
  Eo <- cbind(0, b[["mub_Intercept"]] + b[["mub_x"]] * dd$x, 0.2)
  llo <- sum(ref_dirichlet(dd$Y, sx_softmax(Eo) *
                             summary(fo)$spec_pars["phi", "Estimate"]))
  expect_lt(abs(llo - as.numeric(logLik(fo))), 1e-8 * abs(llo))
})

test_that("categorical() and mixture() refuse a constant on a location", {
  set.seed(20)
  n <- 200
  d <- data.frame(x = stats::rnorm(n))
  d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
  for (dp in c("mub", "muc")) {
    bfo <- do.call(bf, c(list(f ~ x), stats::setNames(list(0), dp)))
    expect_error(frm(bfo, family = categorical(), data = d),
                 "cannot be held at a value", label = dp)
  }
  k <- stats::rbinom(n, 1, 0.4)
  d$y <- ifelse(k == 1, 2 + d$x, -1) + stats::rnorm(n, 0, 0.7)
  mix <- mixture(gaussian(), gaussian())
  expect_error(frm(bf(y ~ x, mu1 = 0), family = mix, data = d),
               "cannot be held at a value")
  expect_error(frm(bf(y ~ x, mu2 = 0), family = mix, data = d),
               "cannot be held at a value")
  # equating two locations keeps brms's own refusal
  expect_error(frm(bf(y ~ x, mu2 = "mu1"), family = mix, data = d),
               "Equating parameters of class 'mu' is not allowed")
})

test_that("a numeric refcat is read as its string, as brms reads it", {
  dd <- sx_data_dir(21)
  dd$Y <- unname(dd$Y)
  f <- frm(bf(Y ~ x), family = dirichlet(refcat = 2), data = dd)
  expect_identical(rownames(fixef(f))[1:2],
                   c("mu1_Intercept", "mu3_Intercept"))
  expect_error(logistic_normal(refcat = c(1, 2)), "with one string")
})

test_that("the routes a composition has no answer on say so", {
  dd <- sx_data_dir(22)
  fit <- frm(bf(Y ~ x), family = dirichlet(), data = dd)
  # every type and new data reach the one refusal, which points to
  # predict() and not to the ordinal routes
  for (ty in c("response", "pearson", "deviance")) {
    expect_error(residuals(fit, type = ty), "Compare predict[(]fit[)]")
  }
  expect_error(residuals(fit, newdata = dd[1:5, ]), "row of proportions")
})

test_that("dharma_residuals() refuses a composition in words", {
  # it died on the matrix response, "(subscript) logical subscript too
  # long", and residuals() used to point there
  skip_if_not_installed("DHARMa")
  dd <- sx_data_dir(22)
  fit <- frm(bf(Y ~ x), family = dirichlet(), data = dd)
  expect_error(dharma_residuals(fit, nsim = 10),
               "a row of proportions has no such function")
})

test_that("conditional_effects(method = \"predict\") says why it is refused", {
  dd <- sx_data_dir(23)
  fd <- frm(bf(Y ~ x), family = dirichlet(), data = dd)
  expect_error(conditional_effects(fd, "x", method = "posterior_predict"),
               "Use method = \"epred\"", fixed = TRUE)
  dl <- sx_data_ln(23)
  fl <- frm(bf(Y ~ x), family = logistic_normal(), data = dl)
  e <- tryCatch(conditional_effects(fl, "x", method = "posterior_predict"),
                error = function(e) conditionMessage(e))
  expect_match(e, "no closed form")
  expect_false(grepl("epred", e, fixed = TRUE))
  # and the predictor it does point to answers
  ce <- conditional_effects(fl, "x", dpar = "mub")
  expect_s3_class(ce, "frmtmb_conditional_effects")
})
