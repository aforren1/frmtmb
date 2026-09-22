# What the first review of items 2.6d and 2.6f found: two blockers,
# five majors and the minors worth fixing
# (dev/reviews/20260918-shapes.md).
#
# Every assertion here was run against the build the review read and
# recorded failing: dev/shapes-log/seen-failing-punch1.txt for this
# file, and dev/shapes-log/interop-lane-seen-failing.txt for the four
# items whose failure is only visible through a third-party package.

punch_fit <- local({
  cache <- new.env(parent = emptyenv())
  function(key) {
    if (!is.null(cache[[key]])) return(cache[[key]])
    set.seed(20260918)
    n <- 150
    dd <- data.frame(x = rnorm(n), z = rnorm(n),
                     g = factor(rep(1:15, each = 10)))
    dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(15, 0, 0.8)[dd$g], 1)
    dd$y2 <- rnorm(n, 0.3 - 0.2 * dd$x, 1)
    dd$ord <- factor(cut(1 + 0.8 * dd$x + rnorm(n),
                         c(-Inf, -0.3, 0.8, Inf), labels = 1:3),
                     ordered = TRUE)
    out <- switch(key,
      mixed = frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd),
      plain = frm(bf(y ~ x) + gaussian(), data = dd),
      ord = frm(bf(ord ~ x) + cumulative(), data = dd),
      ordcs = frm(bf(ord ~ x + cs(z)) + sratio(), data = dd),
      mv = frm(bf(mvbind(y, y2) ~ x) + gaussian(), data = dd))
    attr(out, "punch_data") <- dd
    cache[[key]] <- out
    out
  }
})

punch_data <- function(key) attr(punch_fit(key), "punch_data")

# BLOCKER 1. insight::get_residuals() must stay a vector.

test_that("insight::get_residuals() is one number per observation", {
  skip_if_not_installed("insight")
  fit <- punch_fit("plain")
  r <- insight::get_residuals(fit)
  # residuals() is brms's n x 4 matrix now, and insight's default
  # handed that matrix straight back: mean(), sd() and qqnorm() on the
  # result all changed value in silence
  expect_null(dim(r))
  expect_length(r, stats::nobs(fit))
  expect_s3_class(r, "insight_residuals")
  expect_equal(as.numeric(r),
               unname(stats::residuals(fit)[, "Estimate"]))
  # the whole matrix is still reachable, as it is for a brmsfit
  expect_equal(dim(attr(r, "full")), c(stats::nobs(fit), 4L))
})

# BLOCKER 2. An unseen level's effect is drawn, not held at zero.

test_that("predict(allow_new_levels = TRUE) widens at an unseen level", {
  fit <- punch_fit("mixed")
  nd <- data.frame(x = c(0, 0), g = factor(c("1", "nosuch")))
  set.seed(11)
  pr <- predict(fit, newdata = nd, allow_new_levels = TRUE, ndraws = 4000)
  known <- pr[1L, "Est.Error"]
  fresh <- pr[2L, "Est.Error"]
  # it used to be bit-identical to the known level, because the unseen
  # level's effect was held at zero
  expect_gt(fresh, known * 1.15)
  # and the target is sqrt(sigma^2 + tau^2 + parameter uncertainty),
  # which frm_linpred(se.fit = TRUE) already held
  vc <- VarCorr(fit)
  tau2 <- vc$g$sd[1L, "Estimate"]^2
  sig2 <- sigma(fit)^2
  expect_equal(unname(fresh), sqrt(sig2 + tau2), tolerance = 0.12)
  expect_equal(unname(known), sqrt(sig2), tolerance = 0.12)
})

test_that("predict() takes brms's sample_new_levels = \"gaussian\"", {
  fit <- punch_fit("mixed")
  nd <- data.frame(x = 0, g = factor("nosuch"))
  set.seed(12)
  a <- predict(fit, newdata = nd, allow_new_levels = TRUE, ndraws = 300,
               sample_new_levels = "gaussian")
  set.seed(12)
  b <- predict(fit, newdata = nd, allow_new_levels = TRUE, ndraws = 300)
  expect_identical(a, b)
  expect_error(predict(fit, newdata = nd, allow_new_levels = TRUE,
                       ndraws = 10, sample_new_levels = "uncertainty"),
               "gaussian")
})

# MAJOR 3. An ordinal fit's thresholds are population-level rows.

test_that("fixef(), vcov() and summary()$fixed carry the thresholds", {
  fit <- punch_fit("ord")
  fe <- fixef(fit)
  # brms's own rows on the same model (dev/shapes-p1-brmsord.R)
  expect_equal(rownames(fe), c("Intercept[1]", "Intercept[2]", "x"))
  V <- vcov(fit)
  expect_equal(dim(V), c(3L, 3L))
  expect_equal(rownames(V), rownames(fe))
  expect_equal(rownames(summary(fit)$fixed), rownames(fe))
  # Est.Error is the delta method over the joint covariance: the
  # threshold is a function of tau_raw, not tau_raw itself
  expect_equal(unname(fe[, "Est.Error"]), unname(sqrt(diag(V))))
  # and the second threshold is the first plus a positive increment
  expect_gt(fe["Intercept[2]", "Estimate"], fe["Intercept[1]", "Estimate"])
  # the printed summary shows all three
  txt <- utils::capture.output(print(fit))
  expect_true(any(grepl("Intercept\\[2\\]", txt)))
})

test_that("cs() coefficients are fixef() rows, after the b block", {
  fit <- punch_fit("ordcs")
  expect_equal(rownames(fixef(fit)),
               c("Intercept[1]", "Intercept[2]", "x", "z[1]", "z[2]"))
  expect_equal(dim(vcov(fit)), c(5L, 5L))
})

# MAJOR 4. The same gap was a wrong standard error at an interop seam.

test_that("get_coef()/get_vcov() carry the ordinal thresholds", {
  skip_if_not_installed("marginaleffects")
  fit <- punch_fit("ord")
  cf <- marginaleffects::get_coef(fit)
  # marginaleffects perturbs this vector one entry at a time, so a
  # parameter missing from it contributes nothing to the standard
  # error: avg_slopes() reported 0.00029 where MASS::polr reports
  # 0.02016 on the middle category
  expect_true(all(c("tau_raw_1", "tau_raw_2") %in% names(cf)))
  V <- marginaleffects::get_vcov(fit)
  expect_equal(rownames(V), names(cf))
  expect_equal(nrow(V), length(cf))
  # set_coef() has to write them back, or the perturbation is a no-op
  cf2 <- cf
  cf2[["tau_raw_1"]] <- cf2[["tau_raw_1"]] + 0.25
  f2 <- marginaleffects::set_coef(fit, unname(cf2))
  expect_equal(unname(marginaleffects::get_coef(f2)), unname(cf2))
})

test_that("insight::get_varcov() covers every get_parameters() row", {
  skip_if_not_installed("insight")
  fit <- punch_fit("ord")
  p <- insight::get_parameters(fit)
  V <- insight::get_varcov(fit)
  expect_equal(nrow(V), nrow(p))
  expect_equal(rownames(V), as.character(p$Parameter))
})

# MAJOR 5. coef() and vcov() name the same parameters.

test_that("coef() uses brms's names, so coef()/vcov() pair up", {
  fit <- punch_fit("plain")
  expect_equal(names(coef(fit)), rownames(vcov(fit)))
  expect_true("Intercept" %in% names(coef(fit)))
  # lmtest::coeftest() intersects the two, so a split dropped the
  # intercept row out of a printed significance table with no warning
  mixed <- punch_fit("mixed")
  expect_equal(colnames(coef(mixed)$g), rownames(vcov(mixed)))
  ordf <- punch_fit("ord")
  expect_equal(names(coef(ordf)), rownames(vcov(ordf)))
  # coef() is still fixef() plus ranef(), under the new names
  expect_equal(coef(mixed)$g[["Intercept"]],
               fixef(mixed)["Intercept", "Estimate"] + ranef(mixed)$g[, 1],
               ignore_attr = TRUE)
})

# MAJOR 6. A multivariate predict() answers for every response.

test_that("predict() on a multivariate fit is brms's n x 4 x nresp", {
  fit <- punch_fit("mv")
  set.seed(13)
  pr <- predict(fit, ndraws = 200)
  expect_equal(length(dim(pr)), 3L)
  expect_equal(dim(pr)[2:3], c(4L, 2L))
  expect_equal(dimnames(pr)[[2L]],
               c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_equal(dimnames(pr)[[3L]], c("y", "y2"))
  # one response at a time is still the matrix
  set.seed(13)
  one <- predict(fit, resp = "y2", ndraws = 200)
  expect_equal(dim(one), c(stats::nobs(fit), 4L))
  # the draws stack the same way
  set.seed(13)
  d <- predict(fit, ndraws = 20, summary = FALSE)
  expect_equal(dim(d), c(20L, stats::nobs(fit), 2L))
  expect_equal(dimnames(d)[[3L]], c("y", "y2"))
})

# The minors the review asked for.

test_that("an ordinal fitted() cannot report a probability outside [0, 1]", {
  fit <- punch_fit("ord")
  ft <- fitted(fit)
  expect_equal(length(dim(ft)), 3L)
  q <- ft[, c("Q2.5", "Q97.5"), ]
  expect_gte(min(q), 0)
  expect_lte(max(q), 1)
})

test_that("empty summary blocks take brms's empty shapes", {
  s <- summary(punch_fit("plain"))
  # brms gives NULL for $random and a 0-row frame for $cor_pars
  expect_null(s$random)
  expect_s3_class(s$cor_pars, "data.frame")
  expect_equal(nrow(s$cor_pars), 0L)
  expect_true(NROW(s$spec_pars) > 0L)
  o <- summary(punch_fit("ord"))
  expect_s3_class(o$spec_pars, "data.frame")
  expect_equal(nrow(o$spec_pars), 0L)
})

test_that("brms leaves the row dimnames of its summaries NULL", {
  fit <- punch_fit("plain")
  expect_null(rownames(fitted(fit)))
  expect_null(rownames(stats::residuals(fit)))
  set.seed(14)
  expect_null(rownames(predict(fit, ndraws = 50)))
  expect_null(rownames(fitted(punch_fit("ord"))))
})

test_that("predict(ndraws = ) refuses a fraction rather than truncating", {
  fit <- punch_fit("plain")
  # as.integer(2.5) is 2, and the check ran on the coerced value
  expect_error(predict(fit, ndraws = 2.5), "whole number")
  expect_error(predict(fit, ndraws = 0), "whole number")
})

test_that("a replicate that draws a non-finite dpar is dropped, not fatal", {
  # BLOCKER 2's draw reaches a parameter the fit barely identifies: on
  # the port's mixture fixture the `(1 | patient)` log standard
  # deviation has variance 7.97e6 in vcov(full = TRUE), so one draw in
  # a handful makes the unseen level's variance exp(2730) and mu1 is
  # Inf. The simulator then failed with "NA in probability vector",
  # which names neither the draw nor the parameter. The summary is now
  # taken over the replicates that are there.
  d <- matrix(rnorm(40), nrow = 8, ncol = 5)
  d[3, ] <- NA_real_
  sm <- frmtmb:::brms_summarize_draws(d)
  expect_false(anyNA(sm))
  expect_equal(unname(sm[, "Estimate"]),
               unname(colMeans(d[-3, ])), tolerance = 1e-12)
  expect_equal(unname(sm[, "Est.Error"]),
               unname(apply(d[-3, ], 2L, stats::sd)), tolerance = 1e-12)
})
