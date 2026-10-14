# Pooled post-processing of a frm_multiple() result (lane surface,
# 2026-10-06): the brms_missings vignette's fixef(), summary() and
# conditional_effects() on a brm_multiple() fit, by Rubin's rules.
#
# Seen to fail on 0.68.1 (rellib-r6): fixef() had no method ("no
# applicable method ... frmtmb_multiple"), the pooled table named the
# intercept "(Intercept)", summary() fell to summary.default's list
# layout, and conditional_effects() refused (dev/surface-repros.R).

imps_gauss <- function(m = 3) {
  set.seed(8)
  n <- 80
  x <- stats::rnorm(n)
  y <- stats::rnorm(n, 1 + 0.5 * x, 1)
  x[sample(n, 15)] <- NA
  lapply(seq_len(m), function(i) {
    xi <- x
    xi[is.na(xi)] <- sample(x[!is.na(x)], sum(is.na(xi)), TRUE)
    data.frame(y = y, x = xi)
  })
}

test_that("fixef() pools fixef() of each imputation under brms's names", {
  imps <- imps_gauss()
  fm <- frm_multiple(bf(y ~ x), family = gaussian(), data = imps)
  fx <- fixef(fm)
  expect_equal(rownames(fx), c("Intercept", "x"))
  expect_equal(colnames(fx), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_equal(rownames(fm$pooled), c("Intercept", "x", "sigma_Intercept"))
  # Rubin's rules by hand over the per-imputation fixef()
  per <- lapply(fm$fits, fixef)
  Q <- sapply(per, function(m) m[, "Estimate"])
  U <- sapply(per, function(m) m[, "Est.Error"]^2)
  m <- length(per)
  tv <- rowMeans(U) + (1 + 1 / m) * apply(Q, 1, stats::var)
  expect_equal(unname(fx[, "Estimate"]), unname(rowMeans(Q)))
  expect_equal(unname(fx[, "Est.Error"]), unname(sqrt(tv)))
  # the same numbers as the table frm_multiple() always printed
  expect_equal(unname(fx[, "Estimate"]), fm$pooled$estimate[1:2])
  expect_equal(unname(fx[, "Est.Error"]), fm$pooled$se[1:2])
  # the interval is a t interval on the pooled df
  q <- stats::qt(0.975, fm$pooled$df[1:2])
  expect_equal(unname(fx[, "Q97.5"]),
               fm$pooled$estimate[1:2] + q * fm$pooled$se[1:2])
  expect_equal(fixef(fm, pars = "x"), fx["x", , drop = FALSE])
  expect_error(fixef(fm, summary = FALSE), "summary = FALSE")
  expect_error(fixef(fm, nosuch = 1), "has no argument `nosuch`")
})

test_that("fixef() agrees with mice's pooling of lm()", {
  skip_if_not_installed("mice")
  data("nhanes", package = "mice")
  imp <- mice::mice(nhanes, m = 5, print = FALSE, seed = 1)
  fm <- frm_multiple(bmi ~ age * chl, data = imp)
  ref <- summary(mice::pool(with(imp, stats::lm(bmi ~ age * chl))))
  fx <- fixef(fm)
  expect_equal(rownames(fx), c("Intercept", "age", "chl", "age:chl"))
  # gaussian ML beta-hat is OLS, so the pooled estimates agree to the
  # optimizer's precision; the SEs differ by the ML/REML residual
  # variance convention, as test-multiple-pooling.R says
  expect_equal(unname(fx[, "Estimate"]), ref$estimate, tolerance = 1e-5)
  expect_lt(max(abs(fx[, "Est.Error"] / ref$std.error - 1)), 0.25)
})

test_that("summary() prints the pooled blocks in brms's layout", {
  fm <- frm_multiple(bf(y ~ x), family = gaussian(), data = imps_gauss())
  s <- summary(fm)
  expect_s3_class(s, "summary.frmtmb_multiple")
  expect_equal(rownames(s$fixed), c("Intercept", "x"))
  expect_equal(names(s$fixed), c("Estimate", "Est.Error", "l-95% CI",
                                 "u-95% CI", "df", "fmi"))
  expect_equal(s$fixed$Estimate, unname(fixef(fm)[, "Estimate"]))
  # sigma, pooled on the log scale and reported on its own
  expect_equal(rownames(s$spec_pars), "sigma")
  expect_equal(s$spec_pars$Estimate, exp(fm$pooled["sigma_Intercept",
                                                    "estimate"]))
  expect_lt(s$spec_pars[["l-95% CI"]], s$spec_pars$Estimate)
  expect_output(print(s), "Regression Coefficients:")
  expect_output(print(s), "Rubin's rules")
  expect_output(print(s), "Further Distributional Parameters:")
  expect_error(summary(fm, waic = TRUE), "cannot honor `waic`")
  expect_error(summary(fm, robust = TRUE), "robust = TRUE")
})

test_that("conditional_effects() pools each grid point by Rubin's rules", {
  fm <- frm_multiple(bf(y ~ x), family = gaussian(), data = imps_gauss())
  ce <- conditional_effects(fm, "x")
  expect_s3_class(ce, "frmtmb_conditional_effects")
  d <- ce$x
  # by hand: every imputation's curve on the first one's grid
  nd <- d["x"]
  P <- sapply(fm$fits, function(f) {
    frm_linpred(f, newdata = nd, type = "link", se.fit = TRUE)$fit
  })
  S <- sapply(fm$fits, function(f) {
    frm_linpred(f, newdata = nd, type = "link", se.fit = TRUE)$se.fit
  })
  m <- ncol(P)
  tv <- rowMeans(S^2) + (1 + 1 / m) * apply(P, 1, stats::var)
  expect_equal(d$estimate__, unname(rowMeans(P)))
  expect_equal(d$se__, unname(sqrt(tv)))
  expect_true(all(d$lower__ < d$estimate__ & d$estimate__ < d$upper__))
  # wider than one imputation's band: the between part is in it
  d1 <- conditional_effects(fm$fits[[1]], "x")$x
  expect_true(all(d$se__ >= d1$se__))
  # the bands that do not pool are refused by name
  expect_error(conditional_effects(fm, "x", method = "predict"),
               "cannot honor method = \"predict\"", fixed = TRUE)
  expect_error(conditional_effects(fm, "x", band = "boot"),
               "cannot honor band = \"boot\"", fixed = TRUE)
})

test_that("an ordinal frm_multiple() pools its category probabilities", {
  set.seed(7)
  n <- 150
  ord <- lapply(1:3, function(i) {
    x <- stats::rnorm(n)
    data.frame(x = x, y = cut(x + stats::rlogis(n), c(-Inf, -1, 0, 1, Inf),
                              labels = FALSE))
  })
  fo <- frm_multiple(y ~ x, family = cumulative(), data = ord)
  fx <- fixef(fo)
  expect_equal(rownames(fx), c("Intercept[1]", "Intercept[2]",
                               "Intercept[3]", "x"))
  co <- conditional_effects(fo, "x")
  d <- co[[1L]]
  expect_true(all(d$lower__ >= 0 & d$upper__ <= 1))
  expect_true(all(d$lower__ <= d$estimate__ & d$estimate__ <= d$upper__))
  # each pooled probability is the mean of the imputations'
  # probabilities at that grid point, so the categories sum to 1 as each
  # imputation's do. Lane postfix (2026-10-13): on 0.70.0 (rellib-r8) it
  # was the back-transformed mean logit, and the categories of a
  # cumulative fit over 3 imputations summed to 0.99669 to 0.99983
  # (dev/postfix-rev2-rubin.R)
  per <- sapply(fo$fits, function(f) {
    dj <- conditional_effects(f, "x", int_conditions = list(x = unique(d$x)))
    dj <- dj[[1L]]
    dj$estimate__[match(paste(d$x, d$cats__), paste(dj$x, dj$cats__))]
  })
  expect_equal(d$estimate__, rowMeans(per))
  s <- tapply(d$estimate__, d$x, sum)
  expect_equal(as.vector(s), rep(1, length(s)))
  # fitted() pools by the same rule: on the curve's grid it is the curve,
  # its categories sum to 1, and each estimate is inside its interval
  xs <- sort(unique(d$x))
  fi <- fitted(fo, newdata = data.frame(x = xs))
  i <- match(d$x, xs)
  k <- match(as.character(d$cats__), as.character(1:4))
  expect_equal(fi[cbind(i, 1L, k)], d$estimate__)
  expect_equal(fi[cbind(i, 3L, k)], d$lower__)
  expect_equal(fi[cbind(i, 4L, k)], d$upper__)
  expect_equal(unname(apply(fi[, "Estimate", ], 1, sum)), rep(1, length(xs)))
  expect_true(all(fi[, "Q2.5", ] <= fi[, "Estimate", ] &
                    fi[, "Estimate", ] <= fi[, "Q97.5", ]))
  expect_true(all(fi[, 3:4, ] >= 0 & fi[, 3:4, ] <= 1))
})

test_that("plot() of a frm_multiple() result names what to use instead", {
  fm <- frm_multiple(bf(y ~ x), family = gaussian(), data = imps_gauss())
  expect_error(plot(fm, variable = "^b", regex = TRUE),
               "brms's display is the posterior draws", fixed = TRUE)
})

# Lane postfix (2026-10-13): vcov(), fitted(), coef(), ranef(),
# VarCorr(), emmeans(), nobs(), formula() and family() pool or pass
# through, and the methods with no Rubin's-rules meaning are refused by
# name. Seen to fail on 0.70.0 (rellib-r8, dev/postfix-log/
# multiple-methods-before.txt): coef() and fitted() returned NULL in
# silence (stats' default methods read `$coefficients` and
# `$fitted.values`), vcov(), ranef(), VarCorr(), predict() and logLik()
# had no method, and emmeans() could not handle the class.

imps_pois <- function() {
  set.seed(11)
  n <- 120
  d <- data.frame(x1 = stats::rnorm(n),
                  x2 = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- stats::rpois(n, exp(0.4 + 0.5 * d$x1 +
                               c(a = 0, b = 0.3, c = -0.2)[d$x2]))
  d$x1[sample(n, 30)] <- NA
  d
}

test_that("vcov(), fixef() and fitted() agree with glm() pooled by mice", {
  skip_if_not_installed("mice")
  imp <- mice::mice(imps_pois(), m = 5, print = FALSE, seed = 12)
  fm <- frm_multiple(y ~ x1 + x2, family = poisson(), data = imp)
  mi <- with(imp, stats::glm(y ~ x1 + x2, family = stats::poisson()))
  gl <- mice::getfit(mi)
  pd <- mice::pool(mi)$pooled
  # a poisson glm()'s covariance is the maximum-likelihood one frm()
  # has, so the two pool the same numbers; the gap is the optimizers',
  # judged against the pooled standard error of each coefficient
  fx <- fixef(fm)
  expect_lt(max(abs(fx[, "Estimate"] - pd$estimate) / sqrt(pd$t)), 1e-3)
  expect_lt(max(abs(fx[, "Est.Error"] / sqrt(pd$t) - 1)), 1e-3)
  # vcov() is Rubin's total covariance, off-diagonal entries included
  W <- Reduce(`+`, lapply(gl, stats::vcov)) / 5
  Tv <- W + (1 + 1 / 5) * stats::cov(t(sapply(gl, stats::coef)))
  V <- vcov(fm)
  expect_equal(dimnames(V), list(rownames(fx), rownames(fx)))
  expect_lt(max(abs(V - Tv)) / max(abs(Tv)), 1e-3)
  expect_equal(sqrt(diag(V)), fx[, "Est.Error"])
  expect_equal(vcov(fm, correlation = TRUE), stats::cov2cor(V))
  # fitted() on new rows: Rubin's rules over each imputation's
  # predict.glm() on the link scale, carried back through exp(), and
  # mice's own pooling of the same linear predictor (fix round 1: it
  # pooled on the response scale)
  nd <- data.frame(x1 = c(-1, 0, 1.5),
                   x2 = factor(c("a", "b", "c"), levels = c("a", "b", "c")))
  pg <- lapply(gl, function(g) {
    stats::predict(g, nd, type = "link", se.fit = TRUE)
  })
  qe <- sapply(pg, `[[`, "fit")
  se <- sqrt(rowMeans(sapply(pg, function(p) p$se.fit^2)) +
               (1 + 1 / 5) * apply(qe, 1, stats::var))
  pm <- vapply(seq_len(nrow(nd)), function(i) {
    mice::pool.scalar(qe[i, ], sapply(pg, function(p) p$se.fit[i]^2),
                      n = Inf)$qbar
  }, 0)
  expect_equal(pm, unname(rowMeans(qe)))
  fi <- fitted(fm, newdata = nd)
  expect_equal(dim(fi), c(3L, 4L))
  expect_lt(max(abs(log(fi[, "Estimate"]) - rowMeans(qe)) / se), 1e-3)
  expect_lt(max(abs(fi[, "Est.Error"] / (exp(rowMeans(qe)) * se) - 1)),
            1e-3)
  expect_true(all(fi[, "Q2.5"] > 0))
  # the interval is the link-scale t interval on the Barnard-Rubin df,
  # carried back: the same pooling run by hand on glm()'s predictor
  U <- sapply(pg, function(p) p$se.fit^2)
  rb <- rubin_pool(qe, U, df.residual(fm$fits[[1]]))
  hi <- exp(rb$estimate + stats::qt(0.975, rb$df) * rb$se)
  expect_lt(max(abs(log(fi[, "Q97.5"]) - log(hi)) / rb$se), 1e-3)
  # in sample, each imputation predicts its own rows
  expect_equal(dim(fitted(fm)), c(120L, 4L))
  expect_identical(nobs(fm), 120)
})

# Fix round 1 (review item 4, dev/postfix-rev-rubin.R section 2): on a
# bernoulli model fitted() pooled on the response scale gave
# Q2.5 -0.1668 at x = -1.906 where the pooled curve gave a lower bound
# of 0.0040. Both now pool on the link scale.
test_that("fitted() and the pooled curve agree, inside the support", {
  skip_if_not_installed("mice")
  set.seed(43)
  n <- 40
  x <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, stats::plogis(2 + 2 * x))
  x[sample(n, 10)] <- NA
  imp <- mice::mice(data.frame(y, x), m = 5, seed = 44, print = FALSE)
  fb <- frm_multiple(y ~ x, family = bernoulli(), data = imp)
  ce <- conditional_effects(fb, "x")$x
  fi <- fitted(fb, newdata = data.frame(x = ce$x))
  expect_equal(unname(fi[, "Estimate"]), ce$estimate__)
  expect_equal(unname(fi[, "Q2.5"]), ce$lower__)
  expect_equal(unname(fi[, "Q97.5"]), ce$upper__)
  expect_true(all(fi[, "Q2.5"] > 0 & fi[, "Q97.5"] < 1))
  # the linear predictor pools as it is: no transform
  fl <- fitted(fb, newdata = data.frame(x = ce$x), scale = "linear")
  expect_equal(stats::plogis(unname(fl[, "Estimate"])), ce$estimate__)
})

test_that("emmeans() of a pooled fit is emmeans's own mira pooling", {
  skip_if_not_installed("mice")
  skip_if_not_installed("emmeans")
  imp <- mice::mice(imps_pois(), m = 5, print = FALSE, seed = 12)
  fm <- frm_multiple(y ~ x1 + x2, family = poisson(), data = imp)
  mi <- with(imp, stats::glm(y ~ x1 + x2, family = stats::poisson()))
  e1 <- as.data.frame(emmeans::emmeans(fm, ~ x2))
  e2 <- as.data.frame(emmeans::emmeans(mi, ~ x2))
  expect_lt(max(abs(e1$emmean - e2$emmean) / e2$SE), 1e-3)
  expect_lt(max(abs(e1$SE / e2$SE - 1)), 1e-3)
  # the Barnard-Rubin df of each mean, on emmeans's asymptotic
  # complete-data reference
  expect_lt(max(abs(e1$df / e2$df - 1)), 1e-2)
  c1 <- as.data.frame(pairs(emmeans::emmeans(fm, ~ x2)))
  c2 <- as.data.frame(pairs(emmeans::emmeans(mi, ~ x2)))
  expect_lt(max(abs(c1$estimate - c2$estimate) / c2$SE), 1e-3)
  expect_lt(max(abs(c1$df / c2$df - 1)), 1e-2)
})

imps_mixed <- function() {
  set.seed(21)
  n <- 300
  g <- factor(rep(1:20, length.out = n))
  u0 <- stats::rnorm(20, 0, 0.8)
  u1 <- stats::rnorm(20, 0, 0.4)
  x <- stats::rnorm(n)
  y <- 1 + u0[g] + (0.5 + u1[g]) * x + stats::rnorm(n)
  x[sample(n, 60)] <- NA
  lapply(1:3, function(i) {
    xi <- x
    xi[is.na(xi)] <- sample(x[!is.na(x)], sum(is.na(xi)), TRUE)
    data.frame(y = y, x = xi, g = g)
  })
}

test_that("ranef(), coef() and VarCorr() pool a mixed model", {
  fm <- frm_multiple(y ~ x + (1 + x | g), data = imps_mixed())
  m <- fm$m
  per <- lapply(fm$fits, ranef, condVar = TRUE)
  re <- ranef(fm, condVar = TRUE)
  expect_s3_class(re, "ranef_frmtmb")
  expect_equal(dimnames(re$g), dimnames(per[[1]]$g))
  E <- sapply(per, function(r) r$g[, "x"])
  S <- sapply(per, function(r) attr(r$g, "condSD")[, "x"])
  expect_equal(re$g[, "x"], rowMeans(E))
  expect_equal(attr(re$g, "condSD")[, "x"],
               sqrt(rowMeans(S^2) + (1 + 1 / m) * apply(E, 1, stats::var)))
  # coef() is fixef() plus ranef() in each imputation, so its mean is
  # the pooled fixef() plus the pooled ranef()
  cf <- coef(fm)
  expect_equal(cf$g[["x"]], fixef(fm)["x", "Estimate"] + re$g[, "x"],
               ignore_attr = TRUE)
  # VarCorr() pools on the scales of the pooled table frm_multiple()
  # prints, so the two agree; the intervals differ only by the delta
  # method the two run through, judged against the interval's width
  vc <- VarCorr(fm)
  pv <- fm$pooled_varcorr
  sdv <- pv[pv$type == "sd", ]
  expect_equal(unname(vc$g$sd[, "Estimate"]), sdv$estimate)
  expect_lt(max(abs(vc$g$sd[, "Q2.5"] - sdv$lwr) / (sdv$upr - sdv$lwr)),
            1e-4)
  co <- pv[pv$type == "cor", ]
  expect_equal(vc$g$cor[1, "Estimate", 2], co$estimate)
  expect_identical(unname(vc$g$cor[1, c("Estimate", "Est.Error"), 1]),
                   c(1, 0))
  expect_true(vc$residual__$sd[, "Q2.5"] < vc$residual__$sd[, "Estimate"])
  expect_error(VarCorr(fm, robust = TRUE), "robust = TRUE")
  expect_error(ranef(fm, summary = FALSE), "summary = FALSE")
})

test_that("a level missing from one imputation is refused, not pooled", {
  imps <- imps_mixed()
  imps[[2]] <- imps[[2]][imps[[2]]$g != "20", ]
  imps[[2]]$g <- droplevels(imps[[2]]$g)
  fm <- frm_multiple(y ~ x + (1 | g), data = imps)
  expect_error(ranef(fm), "imputation 2 gives a differently shaped")
  expect_error(coef(fm), "imputation 2 gives a differently shaped")
  expect_error(fitted(fm), "imputation 2 gives a differently shaped")
  expect_error(nobs(fm), "different numbers of observations")
  # the same data with every level present pools
  fm3 <- frm_multiple(y ~ x + (1 | g), data = imps_mixed())
  expect_equal(nrow(ranef(fm3)$g), 20L)
  # a fixed-effect level missing from the first imputation: 0.70.0
  # stopped inside vapply(), "values must be length 4"
  imf <- imps_gauss()
  imf[[1]]$f <- factor(rep(c("a", "b"), 40))
  imf[[2]]$f <- factor(rep(c("a", "b", "c"), length.out = 80))
  imf[[3]]$f <- imf[[2]]$f
  expect_error(frm_multiple(bf(y ~ x + f), family = gaussian(), data = imf),
               "imputation 2 has other coefficients than the first")
})

test_that("formula(), family() pass through; ndraws is said once", {
  fm <- frm_multiple(bf(y ~ x), family = gaussian(), data = imps_gauss())
  expect_equal(formula(fm), formula(fm$fits[[1]]))
  expect_identical(family(fm)$family, "gaussian")
  msgs <- character()
  fi <- withCallingHandlers(fitted(fm, ndraws = 10),
                            message = function(m) {
                              msgs <<- c(msgs, conditionMessage(m))
                              invokeRestart("muffleMessage")
                            })
  expect_identical(fi, fitted(fm))
  expect_length(grep("changes nothing", msgs), 1L)
  expect_error(fitted(fm, summary = FALSE), "summary = FALSE")
})

test_that("the methods with no pooled meaning are refused by name", {
  fm <- frm_multiple(bf(y ~ x), family = gaussian(), data = imps_gauss())
  cases <- list(
    list(quote(predict(fm)), "predict[(][)] has no pooled answer",
         "fitted[(][)] gives the pooled expected response"),
    list(quote(residuals(fm)), "residuals[(][)] has no pooled answer",
         "x[$]fits"),
    list(quote(logLik(fm)), "logLik[(][)] .*has no pooled answer",
         "anova[(]m1, m0[)]"),
    list(quote(AIC(fm)), "logLik[(][)] [(]and so AIC", "anova"),
    list(quote(pp_check(fm)), "pp_check[(][)] has no pooled answer",
         "x[$]fits"),
    list(quote(loo(fm)), "loo[(][)] has no pooled answer", "anova"),
    list(quote(waic(fm)), "waic[(][)] has no pooled answer", "anova"),
    list(quote(update(fm, . ~ 1)), "update[(][)] has no pooled answer",
         "frm_multiple[(][)] again"),
    list(quote(stats::model.frame(fm)),
         "model.frame[(][)] has no pooled answer", "x[$]fits"))
  for (cs in cases) {
    msg <- tryCatch(eval(cs[[1]]), error = conditionMessage)
    expect_match(msg, cs[[2]], info = deparse(cs[[1]]))
    expect_match(msg, cs[[3]], info = deparse(cs[[1]]))
  }
  # each of them answers on one imputation's fit, which is what the
  # refusals name
  f1 <- fm$fits[[1]]
  expect_equal(dim(predict(f1, ndraws = 20)), c(80L, 4L))
  expect_equal(dim(residuals(f1)), c(80L, 4L))
  expect_s3_class(logLik(f1), "logLik")
})
