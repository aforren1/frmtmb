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
  # each pooled probability is the back-transformed mean logit of the
  # imputations' probabilities at that grid point, so it lies between
  # their smallest and largest
  per <- sapply(fo$fits, function(f) {
    dj <- conditional_effects(f, "x", int_conditions = list(x = unique(d$x)))
    dj <- dj[[1L]]
    dj$estimate__[match(paste(d$x, d$cats__), paste(dj$x, dj$cats__))]
  })
  expect_equal(d$estimate__, stats::plogis(rowMeans(stats::qlogis(per))))
  expect_true(all(d$estimate__ >= apply(per, 1, min) &
                    d$estimate__ <= apply(per, 1, max)))
})

test_that("plot() of a frm_multiple() result names what to use instead", {
  fm <- frm_multiple(bf(y ~ x), family = gaussian(), data = imps_gauss())
  expect_error(plot(fm, variable = "^b", regex = TRUE),
               "brms's display is the posterior draws", fixed = TRUE)
})
