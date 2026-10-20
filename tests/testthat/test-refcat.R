# brms's `refcat` for categorical() and multinomial(), the bare
# multinomial(), multinomial()'s fitted(), and the linear-scale
# fitted() of a family with several location predictors. What brms
# 2.23.0 does is recorded in dev/catmulti-findings.md, from
# dev/catmulti-brms-probe.R and dev/catmulti-brms-fits.R.

rc_cat_data <- function(seed, n = 300) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  E <- cbind(0, 0.3 + 0.8 * dd$x, -0.4 - 0.5 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  dd$f <- factor(apply(P, 1, function(p) {
    sample(c("a", "b", "c"), 1, prob = p)
  }), levels = c("a", "b", "c"))
  dd
}

rc_mult_data <- function(seed, n = 200, size = 10) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  E <- cbind(0, 0.3 + 0.8 * dd$x, -0.4 - 0.5 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  Yc <- t(vapply(seq_len(n), function(i) {
    as.numeric(stats::rmultinom(1, size, P[i, ]))
  }, numeric(3)))
  colnames(Yc) <- c("y1", "y2", "y3")
  dd$Yc <- Yc
  dd$N <- size
  dd
}

# the largest relative difference, against the larger magnitude
rc_rel <- function(a, b) {
  max(abs(a - b)) / max(abs(a), abs(b))
}

test_that("categorical(refcat =) names the reference as brms does", {
  dd <- rc_cat_data(1)
  f0 <- frm(bf(f ~ x), family = categorical(), data = dd)
  fb <- frm(bf(f ~ x), family = categorical(refcat = "b"), data = dd)
  # brms 2.23.0: dpars mua, muc under refcat = "b"
  expect_identical(rownames(fixef(fb)),
                   c("mua_Intercept", "muc_Intercept", "mua_x", "muc_x"))
  ll <- abs(as.numeric(logLik(f0)))
  expect_lt(abs(as.numeric(logLik(f0) - logLik(fb))), 1e-8 * ll)
  # a reparameterization: mua = -mub and muc = muc - mub under reference a
  e0 <- fixef(f0)[, "Estimate"]
  eb <- fixef(fb)[, "Estimate"]
  se <- fixef(fb)[, "Est.Error"]
  expect_lt(abs(eb[["mua_x"]] + e0[["mub_x"]]), 1e-3 * se[["mua_x"]])
  expect_lt(abs(eb[["muc_x"]] - (e0[["muc_x"]] - e0[["mub_x"]])),
            1e-3 * se[["muc_x"]])
  # the category probabilities do not depend on the reference
  P0 <- fitted(f0)[, "Estimate", ]
  Pb <- fitted(fb)[, "Estimate", ]
  expect_identical(dimnames(Pb), dimnames(P0))
  expect_lt(rc_rel(P0, Pb), 1e-5)
  # a category's own formula is addressed by its new name
  fz <- frm(bf(f ~ x, mua ~ z), family = categorical(refcat = "b"),
            data = dd)
  expect_identical(names(fixef_by_dpar(fz)), c("mua", "muc"))
  expect_named(fixef_by_dpar(fz)$mua, c("(Intercept)", "z"))
})

test_that("a brms categorical() object carries its refcat", {
  skip_if_not_installed("brms")
  dd <- rc_cat_data(2)
  # 0.71.0 dropped it and fitted reference a under the names mub, muc
  fb <- frm(bf(f ~ x), family = brms::categorical(refcat = "b"), data = dd)
  expect_identical(names(fixef_by_dpar(fb)), c("mua", "muc"))
  fn <- frm(bf(f ~ x), family = brmsfamily("categorical", refcat = "c"),
            data = dd)
  expect_identical(names(fixef_by_dpar(fn)), c("mua", "mub"))
})

test_that("refcat is checked as brms checks it", {
  dd <- rc_cat_data(3)
  expect_error(frm(bf(f ~ x), family = categorical(refcat = "zz"),
                   data = dd),
               "The reference response category must be one of a, b, c")
  expect_error(categorical(refcat = c("a", "b")), "with one string")
  expect_error(multinomial(refcat = c("a", "b")), "with one string")
  # brms's as_one_character(): a number names the category "2"
  dd$fi <- as.integer(dd$f)
  fk <- frm(bf(fi ~ x), family = categorical(K = 3, refcat = 2), data = dd)
  expect_identical(names(fixef_by_dpar(fk)), c("mu1", "mu3"))
  # brms's make_stan_names(): dots and underscores leave a dpar name
  dd$g <- factor(c("p.q", "r_s", "t")[as.integer(dd$f)],
                 levels = c("p.q", "r_s", "t"))
  fg <- frm(bf(g ~ x), family = categorical(refcat = "r_s"), data = dd)
  expect_identical(names(fixef_by_dpar(fg)), c("mupq", "mut"))
})

test_that("refcat = NA predicts every category, and says what is lost", {
  # brms allows it and leans on its priors; brms 2.23.0's own fit of
  # this model with its default priors did not mix (R-hat 1.71, the
  # slopes near 519 with sd 696; dev/catmulti-out/brms-fits.txt). At the
  # maximum the likelihood is flat along a common shift, which the SE
  # check names, and the category probabilities are still determined.
  dd <- rc_cat_data(4)
  se_na <- "Standard errors are not available"
  allow_warnings(
    fn <- frm(bf(f ~ x), family = categorical(refcat = NA), data = dd),
    se_na, require = se_na)
  expect_identical(names(fixef_by_dpar(fn)), c("mua", "mub", "muc"))
  f0 <- frm(bf(f ~ x), family = categorical(), data = dd)
  ll <- abs(as.numeric(logLik(f0)))
  expect_lt(abs(as.numeric(logLik(fn) - logLik(f0))), 1e-6 * ll)
  P0 <- frm_linpred(f0, type = "response")
  Pn <- frm_linpred(fn, type = "response")
  expect_lt(rc_rel(P0, Pn), 1e-4)
  # a prior on each intercept and slope identifies the MAP fit
  pr <- c(set_prior("normal(0, 2)", class = "Intercept", dpar = "mua"),
          set_prior("normal(0, 2)", class = "Intercept", dpar = "mub"),
          set_prior("normal(0, 2)", class = "Intercept", dpar = "muc"),
          set_prior("normal(0, 2)", class = "b", dpar = "mua"),
          set_prior("normal(0, 2)", class = "b", dpar = "mub"),
          set_prior("normal(0, 2)", class = "b", dpar = "muc"))
  fp <- frm(bf(f ~ x), family = categorical(refcat = NA), data = dd,
            prior = pr)
  expect_true(all(is.finite(fixef(fp)[, "Est.Error"])))
  # a symmetric prior centers the shift: the predictors sum to zero
  b <- fixef_by_dpar(fp)
  expect_lt(abs(b$mua[["x"]] + b$mub[["x"]] + b$muc[["x"]]),
            1e-3 * fixef(fp)["mua_x", "Est.Error"])
})

test_that("simulate() draws the reference category's share", {
  dd <- rc_cat_data(5, n = 4000)
  fb <- frm(bf(f ~ 1), family = categorical(refcat = "b"), data = dd)
  p <- frm_linpred(fb, type = "response")[1L, ]
  set.seed(51)
  s <- simulate(fb, nsim = 5)
  tab <- table(factor(unlist(s), levels = names(p)))
  m <- sum(tab)
  # each share within 5 binomial standard errors of its probability
  expect_true(all(abs(tab / m - p) < 5 * sqrt(p * (1 - p) / m)))
})

test_that("multinomial() reads its categories from the response", {
  dd <- rc_mult_data(6)
  fk <- frm(bf(Yc | trials(N) ~ x), family = multinomial(K = 3), data = dd)
  fm <- frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dd)
  # brms 2.23.0 names them by the columns: muy2, muy3
  expect_identical(names(fixef_by_dpar(fm)), c("muy2", "muy3"))
  expect_identical(names(fixef_by_dpar(fk)), c("mu2", "mu3"))
  ll <- abs(as.numeric(logLik(fk)))
  expect_lt(abs(as.numeric(logLik(fk) - logLik(fm))), 1e-10 * ll)
  # cbind() names the categories too, and an unnamed matrix numbers them
  fc <- frm(bf(cbind(y1, y2, y3) | trials(N) ~ x),
            family = multinomial(),
            data = cbind(dd, as.data.frame(unclass(dd$Yc))))
  expect_identical(names(fixef_by_dpar(fc)), c("muy2", "muy3"))
  du <- dd
  du$Yc <- unname(dd$Yc)
  fu <- frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = du)
  expect_identical(names(fixef_by_dpar(fu)), c("mu2", "mu3"))
  skip_if_not_installed("brms")
  fb <- frm(bf(Yc | trials(N) ~ x), family = brms::multinomial(), data = dd)
  expect_identical(names(fixef_by_dpar(fb)), c("muy2", "muy3"))
})

test_that("multinomial(refcat =) reparameterizes, as in brms", {
  dd <- rc_mult_data(7)
  f0 <- frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dd)
  f2 <- frm(bf(Yc | trials(N) ~ x), family = multinomial(refcat = "y2"),
            data = dd)
  expect_identical(names(fixef_by_dpar(f2)), c("muy1", "muy3"))
  ll <- abs(as.numeric(logLik(f0)))
  expect_lt(abs(as.numeric(logLik(f0) - logLik(f2))), 1e-8 * ll)
  e0 <- fixef(f0)[, "Estimate"]
  e2 <- fixef(f2)[, "Estimate"]
  expect_lt(abs(e2[["muy1_x"]] + e0[["muy2_x"]]),
            1e-3 * fixef(f2)["muy1_x", "Est.Error"])
  # refcat = NA leaves the common shift to the SE check, as categorical()
  se_na <- "Standard errors are not available"
  allow_warnings(
    fn <- frm(bf(Yc | trials(N) ~ x), family = multinomial(refcat = NA),
              data = dd),
    se_na, require = se_na)
  expect_identical(names(fixef_by_dpar(fn)), c("muy1", "muy2", "muy3"))
  expect_lt(abs(as.numeric(logLik(f0) - logLik(fn))), 1e-6 * ll)
  # simulated counts keep the reference's share and the column names
  set.seed(71)
  s <- simulate(f2, nsim = 1)[[1L]]
  expect_identical(colnames(s), c("y1", "y2", "y3"))
  expect_true(all(rowSums(s) == 10))
})

test_that("multinomial() refusals name the cause", {
  dd <- rc_mult_data(8)
  # the link comes first, as in brms
  expect_error(multinomial(3), "multinomial[(]K = 3[)]")
  dd$Y1 <- matrix(dd$Yc[, 1], ncol = 1)
  expect_error(frm(bf(Y1 | trials(N) ~ x), family = multinomial(),
                   data = dd),
               "At least 2 response categories are required")
  expect_error(frm(bf(Yc | trials(N) ~ x), family = multinomial(K = 4),
                   data = dd),
               "n x 4")
  expect_error(multinomial(K = 1), "at least two categories")
})

test_that("fitted() on multinomial() is brms's n x 4 x K of counts", {
  dd <- rc_mult_data(9)
  # 0.71.0 refused it ("declares no mean") on the K = spelling too
  fk <- frm(bf(Yc | trials(N) ~ x), family = multinomial(K = 3), data = dd)
  expect_identical(dim(fitted(fk)), c(200L, 4L, 3L))
  fm <- frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dd)
  fv <- fitted(fm)
  # brms 2.23.0: 200 x 4 x 3, layers P(Y = y1) ..., rows unnamed
  expect_identical(dim(fv), c(200L, 4L, 3L))
  expect_identical(dimnames(fv),
                   list(NULL, c("Estimate", "Est.Error", "Q2.5", "Q97.5"),
                        c("P(Y = y1)", "P(Y = y2)", "P(Y = y3)")))
  # the expected counts: probability times trials, summing to trials
  est <- fv[, "Estimate", ]
  expect_lt(max(abs(rowSums(est) - dd$N)) / max(dd$N), 1e-12)
  b <- fixef_by_dpar(fm)
  E <- cbind(0, b$muy2[[1]] + b$muy2[[2]] * dd$x,
             b$muy3[[1]] + b$muy3[[2]] * dd$x)
  P <- exp(E) / rowSums(exp(E))
  expect_lt(rc_rel(est, P * dd$N), 1e-10)
  # Est.Error is the delta method through the softmax, against vcov()
  cn <- c("muy2_Intercept", "muy2_x", "muy3_Intercept", "muy3_x")
  V <- vcov(fm)[cn, cn]
  i <- 17L
  xi <- c(1, dd$x[i])
  p <- P[i, ]
  J <- matrix(0, 3, 4)
  # dP_k / d eta_j = p_k (1{k = j} - p_j), j over categories 2 and 3
  for (j in 2:3) {
    d <- p * ((seq_len(3) == j) - p[j])
    J[, (j - 2) * 2 + 1:2] <- outer(d, xi)
  }
  se_ref <- sqrt(diag(J %*% V %*% t(J))) * dd$N[i]
  expect_lt(rc_rel(fv[i, "Est.Error", ], se_ref), 1e-4)
  # a Wald bound stays between 0 and the row's trials
  expect_true(all(fv[, c("Q2.5", "Q97.5"), ] >= 0))
  expect_true(all(sweep(fv[, "Q97.5", ], 1L, dd$N, "<=")))
  # newdata reads its own trials, as brms does
  nd <- data.frame(x = c(-1, 0, 1), N = c(5L, 10L, 20L))
  fn <- fitted(fm, newdata = nd)
  expect_lt(max(abs(rowSums(fn[, "Estimate", ]) - nd$N)) / 20, 1e-12)
})

test_that("fitted(scale = 'linear') stacks every location predictor", {
  dd <- rc_cat_data(10)
  # the default reference, so that 0.71.0 fails on the shape, not the
  # argument: it returned the first predictor alone, 300 x 4
  fb <- frm(bf(f ~ x), family = categorical(), data = dd)
  fl <- fitted(fb, scale = "linear")
  # brms 2.23.0: 300 x 4 x 2, layers eta1, eta2
  expect_identical(dim(fl), c(300L, 4L, 2L))
  expect_identical(dimnames(fl)[[3L]], c("eta1", "eta2"))
  for (k in 1:2) {
    lp <- frm_linpred(fb, dpar = c("mub", "muc")[k], se.fit = TRUE)
    expect_lt(rc_rel(fl[, "Estimate", k], lp$fit), 1e-12)
    expect_lt(rc_rel(fl[, "Est.Error", k], lp$se.fit), 1e-4)
  }
  # without a reference category there are K layers
  na_w <- c("Standard errors are not available",
            "a direction the fit does not determine")
  allow_warnings({
    fn <- frm(bf(f ~ x), family = categorical(refcat = NA), data = dd)
    fnl <- fitted(fn, scale = "linear")
  }, na_w)
  expect_identical(dim(fnl), c(300L, 4L, 3L))
  # a mixture's component predictors are stacked unnamed, as in brms
  set.seed(101)
  dm <- data.frame(x = stats::rnorm(300))
  zz <- stats::rbinom(300, 1, 0.4)
  dm$y <- ifelse(zz == 1, -2 + dm$x, 2 - 0.5 * dm$x) +
    stats::rnorm(300, 0, 0.7)
  fx <- frm(bf(y ~ x), family = mixture(gaussian, gaussian), data = dm)
  fxl <- fitted(fx, scale = "linear")
  expect_identical(dim(fxl), c(300L, 4L, 2L))
  expect_null(dimnames(fxl)[[3L]])
  # one predictor stays a matrix: a two-category categorical
  d2 <- dd
  d2$f <- factor(ifelse(dd$f == "a", "a", "b"))
  f2 <- suppressMessages(frm(bf(f ~ x), family = categorical(), data = d2))
  expect_identical(dim(fitted(f2, scale = "linear")), c(300L, 4L))
})

test_that("multinomial()'s mean reaches the routes that read one", {
  dd <- rc_mult_data(11)
  fm <- frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dd)
  # brms 2.23.0 draws it per category, trials held at 1, and 300 rows
  # over its 100-point grid; frmtmb's per-category display is the
  # default, as for the simplex families (dev/catmulti-ce-compare.R
  # matches brms's estimates to 0.085 of brms's se)
  expect_message(ce <- conditional_effects(fm), "holding the trials")
  ce <- ce[[1L]]
  expect_identical(levels(ce$cats__), c("y1", "y2", "y3"))
  expect_identical(nrow(ce), 300L)
  b <- fixef_by_dpar(fm)
  x0 <- ce$x[1L]
  e <- c(0, b$muy2[[1]] + b$muy2[[2]] * x0, b$muy3[[1]] + b$muy3[[2]] * x0)
  p <- exp(e) / sum(exp(e))
  got <- ce$estimate__[ce$x == x0]
  expect_lt(rc_rel(got, p), 1e-10)
  expect_error(suppressMessages(conditional_effects(fm, method = "predict")),
               "not available for family 'multinomial'")
  # 0.71.0 refused on the missing mean; with one, the rank code failed
  expect_error(dharma_residuals(fm), "has no meaning for the 'multinomial'")
})

# The review's data (dev/catmulti-rev-data.R, rev_cat_data()), copied so
# its numbers reproduce: K categories with a grouping factor g.
rc_rev_data <- function(seed, K, n = 400, ng = 25) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(sample(paste0("g", seq_len(ng)), n, TRUE)))
  L <- t(chol(0.5 * diag(K - 1L) + 0.3))
  u <- matrix(stats::rnorm(ng * (K - 1L)), ng) %*% t(L)
  gi <- as.integer(dd$g)
  b0 <- seq(0.3, -0.4, length.out = K - 1L)
  b1 <- seq(0.8, -0.6, length.out = K - 1L)
  E <- cbind(0, sapply(seq_len(K - 1L), function(k) {
    b0[k] + b1[k] * dd$x + u[gi, k]
  }))
  P <- exp(E) / rowSums(exp(E))
  dd$y <- factor(apply(P, 1, function(p) sample(letters[1:K], 1, prob = p)),
                 levels = letters[1:K])
  dd
}

test_that("refcat = NA counts the reference model's degrees of freedom", {
  # the lane's first build counted the flat shift: df 6 against 4, AIC
  # 626.69 against 622.69, anova() of two copies of one model Chi Df 2,
  # drop1() Df 3 against 2 (dev/reviews/2026-10-14-catmulti.md, 2)
  se_na <- "Standard errors are not available"
  dd <- rc_rev_data(23, 3, n = 300)
  allow_warnings(
    fn <- frm(bf(y ~ x), family = categorical(refcat = NA), data = dd),
    se_na, require = "common shift of every category's predictor")
  fa <- frm(bf(y ~ x), family = categorical(), data = dd)
  expect_identical(attr(logLik(fa), "df"), 4L)
  expect_identical(attr(logLik(fn), "df"), 4L)
  expect_lt(abs(AIC(fn) - AIC(fa)), 1e-8 * AIC(fa))
  expect_lt(abs(AIC(fa) - 622.6884) / 622.6884, 1e-6)
  expect_identical(df.residual(fn), df.residual(fa))
  av <- anova(fa, fn)
  expect_equal(av[["Chi Df"]][2L], 0)
  d3 <- rc_rev_data(21, 3)
  d3$z <- stats::rnorm(nrow(d3))
  allow_warnings(
    fn3 <- frm(bf(y ~ x + z), family = categorical(refcat = NA), data = d3),
    se_na)
  fa3 <- frm(bf(y ~ x + z), family = categorical(), data = d3)
  tn <- allow_warnings(drop1(fn3, test = "Chisq"), se_na)
  ta <- drop1(fa3, test = "Chisq")
  expect_identical(tn$Df, ta$Df)
  expect_identical(ta$Df[-1L], c(2, 2))
  expect_lt(max(abs(tn$AIC - ta$AIC) / ta$AIC), 1e-6)
  # a category with its own formula frees the shift only in part: no
  # df is claimed, and the comparisons say why
  allow_warnings(
    fo <- frm(bf(y ~ x, mua ~ 1), family = categorical(refcat = NA),
              data = dd),
    se_na)
  expect_true(is.na(attr(logLik(fo), "df")))
  expect_error(drop1(fo), "only partly free")
  expect_error(anova(fa, fo), "only partly free")
  # dirichlet(refcat = NA) has the same flat shift
  set.seed(24)
  G <- matrix(stats::rgamma(900, shape = 3), 300, 3)
  dd$Y <- G / rowSums(G)
  colnames(dd$Y) <- c("a", "b", "c")
  allow_warnings(
    dn <- frm(bf(Y ~ x), family = dirichlet(refcat = NA), data = dd), se_na)
  d0 <- frm(bf(Y ~ x), family = dirichlet(), data = dd)
  expect_identical(attr(logLik(dn), "df"), attr(logLik(d0), "df"))
})

test_that("a family's print shows its reference category", {
  out <- utils::capture.output(print(multinomial(refcat = "y2")))
  expect_true(any(grepl("Reference category: y2", out, fixed = TRUE)))
  out <- utils::capture.output(print(categorical(refcat = NA)))
  expect_true(any(grepl("Reference category: none (refcat = NA)", out,
                        fixed = TRUE)))
  out <- utils::capture.output(print(categorical(levels = c("a", "b"))))
  expect_true(any(grepl("Reference category: a", out, fixed = TRUE)))
})

test_that("conditional_effects(categorical = TRUE) reads fitted()", {
  # 0.71.0 refused it, "no analytic standard error ... Use band = boot";
  # dev/catmulti-ce-compare.R matches brms 2.23.0's on the same data
  dd <- rc_cat_data(12)
  fb <- frm(bf(f ~ x), family = categorical(refcat = "b"), data = dd)
  ce <- conditional_effects(fb, "x", categorical = TRUE)[[1L]]
  expect_identical(levels(ce$cats__), c("a", "b", "c"))
  nd <- data.frame(x = ce$x[ce$cats__ == "a"])
  fv <- fitted(fb, newdata = nd)
  for (k in c("a", "b", "c")) {
    lay <- paste0("P(Y = ", k, ")")
    expect_lt(rc_rel(ce$estimate__[ce$cats__ == k], fv[, "Estimate", lay]),
              1e-10)
    expect_lt(rc_rel(ce$se__[ce$cats__ == k], fv[, "Est.Error", lay]),
              1e-10)
  }
})

test_that("refcat = NA counts the df of transformed terms too", {
  # round 1 rebuilt the reference model from the MODEL frame, whose
  # columns are y and log(x2), so any term that reads a raw variable
  # lost its df and AIC (dev/reviews/2026-10-14-catmulti.md, 9)
  se_na <- "Standard errors are not available"
  dd <- rc_rev_data(5, 3, n = 400)
  dd$x2 <- exp(stats::rnorm(400))
  forms <- list(y ~ poly(x, 2) + log(x2), y ~ scale(x) + I(x2^2),
                y ~ log(x2))
  for (fm in forms) {
    allow_warnings(
      fn <- frm(bf(fm), family = categorical(refcat = NA), data = dd),
      se_na)
    fa <- frm(bf(fm), family = categorical(), data = dd)
    expect_identical(attr(logLik(fn), "df"), attr(logLik(fa), "df"),
                     label = deparse1(fm))
    expect_true(is.finite(AIC(fn)), label = deparse1(fm))
    expect_lt(abs(AIC(fn) - AIC(fa)), 1e-6 * AIC(fa))
  }
  # a term update() adds is counted on the refit's own data
  allow_warnings({
    f1 <- frm(bf(y ~ x), family = categorical(refcat = NA), data = dd)
    f2 <- update(f1, . ~ . + log(x2), newdata = dd)
  }, se_na)
  fa2 <- frm(bf(y ~ x + log(x2)), family = categorical(), data = dd)
  expect_identical(attr(logLik(f2), "df"), attr(logLik(fa2), "df"))
  expect_equal(anova(fa2, f2)[["Chi Df"]][2L], 0)
})

test_that("the refcat = NA sentence comes only with the shift", {
  # review nit: it attached to any lost SE of an NA fit, an sd at 0 too
  se_na <- "Standard errors are not available"
  dd <- rc_rev_data(6, 3, n = 300)
  allow_warnings(
    fn <- frm(bf(y ~ x), family = categorical(refcat = NA), data = dd),
    se_na)
  shift <- "common shift of every category's predictor"
  all_lost <- stats::setNames(rep("flat", 2),
                              c("mua_(Intercept)", "mub_(Intercept)"))
  expect_false(grepl(shift, se_lost_message(fn, all_lost), fixed = TRUE))
  all_lost <- c(all_lost, "muc_(Intercept)" = "flat")
  expect_true(grepl(shift, se_lost_message(fn, all_lost), fixed = TRUE))
})
