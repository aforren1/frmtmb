# brms's names, orders and layouts on a fit: variables(), ranef(),
# coef(), the one-response draws matrices, summary()'s Formula and Data
# lines, and confint_varcorr()'s block of a || term.
#
# Every expected vector below is brms 2.23.0's own output on the same
# model, recorded by dev/names-brms.R (logs dev/names-out/brms-*.txt;
# dev/names-findings.md), with the names a fit has no counterpart of
# (r_, the centered Intercept, lp__) left out.

layout_data <- function(n = 120, seed = 20261111) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(rep(sprintf("g%02d", 1:12), each = n / 12)),
                  h = factor(rep(1:8, length.out = n)),
                  ga = factor(rep(1:6, length.out = n)),
                  f = factor(rep(c("fa", "fb"), length.out = n)))
  ug <- stats::rnorm(12, 0, 0.5)[d$g]
  d$y <- 1 + 0.5 * d$x + ug + stats::rnorm(n)
  d$ydist <- 1 + 0.5 * d$x + ug + stats::rnorm(n, 0, exp(0.2 * d$z))
  d$y2 <- 0.5 * d$y + stats::rnorm(n)
  d$k <- sample(1:4, n, replace = TRUE)
  d$ymo <- 1 + 0.4 * d$k + stats::rnorm(n)
  d$sz <- stats::runif(n, 0.2, 0.4)
  d$xg <- stats::rnorm(12)[d$g]
  d$wg <- stats::rnorm(12)[d$g]
  d$sg <- 0.25
  d$yme <- 1 + 0.5 * d$xg + 0.3 * d$z - 0.2 * d$wg + stats::rnorm(n)
  d
}

quiet_fit <- function(expr) suppressMessages(allow_at_max(expr))

test_that("variables() of a distributional fit is in brms's order", {
  d <- layout_data()
  fit <- quiet_fit(frm(bf(ydist ~ x + (1 | g), sigma ~ z + (1 | g)),
                       family = gaussian(), data = d))
  # was b_Intercept b_x b_sigma_Intercept b_sigma_z ...
  expect_identical(variables(fit),
                   c("b_Intercept", "b_sigma_Intercept", "b_x", "b_sigma_z",
                     "sd_g__Intercept", "sd_g__sigma_Intercept"))
  pop <- quiet_fit(frm(bf(ydist ~ x, sigma ~ 1), family = gaussian(),
                       data = d))
  expect_identical(variables(pop),
                   c("b_Intercept", "b_sigma_Intercept", "b_x"))
})

test_that("group ids follow the order the dpar formulas are written in", {
  d <- layout_data()
  d$ys <- d$y + stats::rnorm(nrow(d), 0, 0.3)
  # brms numbers the ids h, alpha g, sigma ga for this writing, and
  # sigma ga, alpha g for the other, while its b_ follow the family
  fa <- quiet_fit(
    frm(bf(ys ~ (1 | h), alpha ~ (1 | g), sigma ~ (1 | ga)),
        family = skew_normal(), data = d))
  expect_identical(variables(fa),
                   c("b_Intercept", "b_sigma_Intercept", "b_alpha_Intercept",
                     "sd_h__Intercept", "sd_g__alpha_Intercept",
                     "sd_ga__sigma_Intercept"))
  expect_identical(names(ranef(fa)), c("h", "g", "ga"))
  expect_identical(names(VarCorr(fa)), c("h", "g", "ga"))
  fs <- quiet_fit(
    frm(bf(ys ~ (1 | h), sigma ~ (1 | ga), alpha ~ (1 | g)),
        family = skew_normal(), data = d))
  expect_identical(variables(fs)[4:6],
                   c("sd_h__Intercept", "sd_ga__sigma_Intercept",
                     "sd_g__alpha_Intercept"))
})

test_that("a multivariate fit lists every intercept first, as brms", {
  d <- layout_data()
  fit <- quiet_fit(frm(bf(y ~ x + (1 | g)) + bf(y2 ~ z + (1 | g)),
                       family = gaussian(), data = d))
  expect_identical(variables(fit),
                   c("b_y_Intercept", "b_y2_Intercept", "b_y_x", "b_y2_z",
                     "sd_g__y_Intercept", "sd_g__y2_Intercept", "sigma_y",
                     "sigma_y2"))
  # one coef() and one ranef() entry for g, over both responses
  expect_named(ranef(fit), "g")
  expect_identical(colnames(ranef(fit)$g), c("y_Intercept", "y2_Intercept"))
  expect_named(coef(fit), "g")
  expect_identical(colnames(coef(fit)$g),
                   c("y_Intercept", "y2_Intercept", "y_x", "y2_z"))
})

test_that("me() hyperparameters are listed group by group", {
  d <- layout_data()
  fit <- quiet_fit(
    frm(bf(yme ~ me(xg, sg, gr = g) + me(z, sz) + me(wg, sg, gr = g)),
        family = gaussian(), data = d))
  v <- variables(fit)
  expect_identical(v[startsWith(v, "meanme_")],
                   c("meanme_mexg", "meanme_mewg", "meanme_mez"))
  expect_identical(v[startsWith(v, "sdme_")],
                   c("sdme_mexg", "sdme_mewg", "sdme_mez"))
})

test_that("variables() of a mo() fit lists its simplex as brms does", {
  d <- layout_data()
  fit <- quiet_fit(frm(bf(ymo ~ mo(k)), family = gaussian(), data = d))
  expect_identical(variables(fit),
                   c("b_Intercept", "bsp_mok", "sigma", "simo_mok1[1]",
                     "simo_mok1[2]", "simo_mok1[3]"))
  # the weights are the summary's, and hypothesis() reads them
  w <- summary(fit)$mo[, "Estimate"]
  # the sum of the weights, against 1, so the tolerance is relative
  h <- hypothesis(fit, "simo_mok1[1] + simo_mok1[2] + simo_mok1[3] = 0",
                  class = NULL)
  expect_equal(h$hypothesis$Estimate, 1, tolerance = 1e-12)
  hv <- hypothesis(fit, "simo_mok1[2] = 0", class = NULL)
  # relative to the weight, which all.equal() divides by
  expect_equal(hv$hypothesis$Estimate, unname(w[2]), tolerance = 1e-10)
})

test_that("ranef() and coef() have one entry per grouping factor", {
  d <- layout_data()
  fit <- quiet_fit(frm(bf(ydist ~ x + (1 | g), sigma ~ z + (1 | g)),
                       family = gaussian(), data = d))
  re <- ranef(fit)
  expect_named(re, "g")
  expect_identical(colnames(re$g), c("Intercept", "sigma_Intercept"))
  expect_identical(attr(re$g, "term"), c("1 | g", "sigma: 1 | g"))
  # a block label still reaches its own columns
  expect_identical(colnames(re[["sigma: 1 | g"]]), "sigma_Intercept")
  expect_identical(unname(re[["sigma: 1 | g"]][, 1]),
                   unname(re$g[, "sigma_Intercept"]))
  df <- as.data.frame(re)
  expect_identical(unique(df$grp), c("1 | g", "sigma: 1 | g"))
  co <- coef(fit)
  expect_named(co, "g")
  expect_identical(colnames(co$g),
                   c("Intercept", "sigma_Intercept", "x", "sigma_z"))
  fe <- fixef(fit)[, "Estimate"]
  expect_equal(co$g$sigma_z, rep(unname(fe["sigma_z"]), 12),
               tolerance = 1e-12)
  expect_equal(co$g$sigma_Intercept,
               unname(fe["sigma_Intercept"]) + unname(re$g[, 2]),
               tolerance = 1e-12)

  # two bars on one factor, brms's one g entry
  f2 <- quiet_fit(frm(bf(y ~ x + (1 | g) + (0 + x | g)), family = gaussian(),
                      data = d))
  re2 <- ranef(f2)
  expect_named(re2, "g")
  expect_identical(colnames(re2$g), c("Intercept", "x"))
  expect_identical(colnames(coef(f2)$g), c("Intercept", "x"))
})

test_that("a categorical fit's ranef() and coef() are brms's", {
  d <- layout_data()
  ug <- stats::rnorm(12, 0, 0.5)[d$g]
  eta <- cbind(0, 0.3 + 0.6 * d$x + ug, -0.2 - 0.4 * d$x)
  p <- exp(eta) / rowSums(exp(eta))
  d$yc <- factor(apply(p, 1, function(pp) {
    sample(c("a", "b", "c"), 1, prob = pp)
  }))
  fit <- quiet_fit(frm(bf(yc ~ x + (1 | g)), family = categorical(),
                       data = d))
  expect_identical(variables(fit),
                   c("b_mub_Intercept", "b_muc_Intercept", "b_mub_x",
                     "b_muc_x", "sd_g__mub_Intercept",
                     "sd_g__muc_Intercept"))
  expect_named(ranef(fit), "g")
  expect_identical(colnames(ranef(fit)$g), c("mub_Intercept", "muc_Intercept"))
  expect_named(coef(fit), "g")
  expect_identical(colnames(coef(fit)$g),
                   c("mub_Intercept", "muc_Intercept", "mub_x", "muc_x"))
})

test_that("ranef() of a model without group-level terms stops, as brms", {
  d <- layout_data()
  fit <- quiet_fit(frm(bf(y ~ s(x)), family = gaussian(), data = d))
  expect_error(ranef(fit), "does not contain group-level effects")
  fit0 <- quiet_fit(frm(bf(y ~ x), family = gaussian(), data = d))
  expect_error(ranef(fit0), "does not contain group-level effects")
  # a model WITH group-level terms keeps answering, also when `groups`
  # selects none of them
  f1 <- quiet_fit(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d))
  expect_length(ranef(f1, groups = "nothing"), 0L)
})

test_that("predict(summary = FALSE) carries brms's dimnames", {
  d <- layout_data()
  fit <- quiet_fit(frm(bf(y ~ x), family = gaussian(), data = d))
  set.seed(1)
  p <- predict(fit, summary = FALSE, ndraws = 5)
  expect_identical(dimnames(p), list(NULL, NULL))
})

test_that("summary() prints every dpar formula on its Formula lines", {
  d <- layout_data()
  fit <- quiet_fit(frm(bf(ydist ~ x + (1 | g), sigma ~ z + (1 | g)),
                       family = gaussian(), data = d))
  out <- utils::capture.output(print(summary(fit)))
  i <- grep("^Formula:", out)
  expect_identical(out[i], "Formula: ydist ~ x + (1 | g) ")
  expect_identical(out[i + 1L], "         sigma ~ z + (1 | g)")
  expect_match(out[i + 2L], "^   Data: d [(]Number of observations")
})

test_that("a data frame passed by value is named by a descriptor", {
  d <- layout_data()
  fit <- quiet_fit(do.call(frm, list(formula = bf(y ~ x), data = d,
                                     family = gaussian())))
  out <- utils::capture.output(print(summary(fit)))
  expect_true(any(grepl(paste0("^   Data: <data[.]frame: ", nrow(d), " x ",
                               ncol(d), "> [(]Number"), out)))
  # names and calls keep brms's deparse
  f2 <- quiet_fit(frm(bf(y ~ x), family = gaussian(), data = d[1:60, ]))
  expect_identical(summary(f2)$data_name, "d[1:60, ]")
  # update(newdata = ) records its data the same two ways
  u1 <- quiet_fit(update(f2, newdata = d))
  expect_identical(summary(u1)$data_name, "d")
  u2 <- quiet_fit(do.call(stats::update, list(f2, newdata = d[1:90, ])))
  expect_identical(summary(u2)$data_name,
                   paste0("<data.frame: 90 x ", ncol(d), ">"))
})

test_that("confint_varcorr() names a || block as written", {
  d <- layout_data()
  d$yf <- 1 + 0.5 * d$x + stats::rnorm(12)[d$g] * (1 + (d$f == "fb")) +
    stats::rnorm(nrow(d))
  fit <- quiet_fit(frm(bf(yf ~ x + (1 + f || g)), family = gaussian(),
                       data = d))
  cv <- confint_varcorr(fit)
  expect_identical(unique(cv$block), "1 + f || g")
  # a correlated block keeps its label
  f2 <- quiet_fit(frm(bf(yf ~ x + (1 + x | g)), family = gaussian(),
                      data = d))
  expect_identical(unique(confint_varcorr(f2)$block), "1 + x | g")
  # the pooled table of several imputations names it the same way, and
  # so does its summary
  fm <- quiet_fit(frm_multiple(bf(yf ~ x + (1 + f || g)), family = gaussian(),
                               data = list(d, d)))
  expect_identical(unique(fm$pooled_varcorr$grp), "1 + f || g")
  out <- utils::capture.output(print(summary(fm)))
  expect_true(any(grepl("1 + f || g", out, fixed = TRUE)))
  expect_false(any(grepl("1 + f | g", out, fixed = TRUE)))
  # a gr(cor = FALSE) term as written too
  f3 <- quiet_fit(frm(bf(yf ~ x + (1 + x | gr(g, cor = FALSE))),
                      family = gaussian(), data = d))
  expect_identical(unique(confint_varcorr(f3)$block),
                   "1 + x | gr(g, cor = FALSE)")
})

test_that("insight reads an order = 'mu' ordinal mixture on fixef()'s scale", {
  skip_if_not_installed("insight")
  set.seed(5)
  n <- 400
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  cls <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.5 * d$z))
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  fit <- quiet_fit(frm(bf(y ~ x),
                       family = mixture(cumulative(),
                                        sratio(threshold = "sum_to_zero"),
                                        order = "mu"), data = d))
  p <- insight::get_parameters(fit)
  fx <- fixef(fit)
  thr <- c(paste0("mu1_Intercept[", 1:3, "]"),
           paste0("mu2_Intercept[", 1:3, "]"))
  # both components' thresholds, where the internal tau_raw used to stand
  expect_true(all(thr %in% p$Parameter))
  expect_false(any(grepl("tau_raw", p$Parameter)))
  est <- stats::setNames(p$Estimate, p$Parameter)[thr]
  expect_equal(unname(est), unname(fx[thr, "Estimate"]),
               tolerance = 1e-10)
  V <- insight::get_varcov(fit)
  expect_identical(rownames(V), p$Parameter)
  expect_equal(sqrt(diag(V))[thr], fx[thr, "Est.Error"],
               tolerance = 1e-6,
               ignore_attr = TRUE)
})

test_that("categorical() drops an unused level and K takes the data's names", {
  set.seed(3164)
  n <- 90
  d <- data.frame(x = stats::rnorm(n))
  d$yc <- factor(sample(c("a", "b", "c"), n, TRUE),
                 levels = c("a", "b", "c", "unused"))
  # brms 2.23.0: ncat 3, b_mub_Intercept b_muc_Intercept b_mub_x b_muc_x
  fit <- quiet_fit(frm(bf(yc ~ x), family = categorical(), data = d))
  expect_identical(variables(fit),
                   c("b_mub_Intercept", "b_muc_Intercept", "b_mub_x",
                     "b_muc_x"))
  fk <- quiet_fit(frm(bf(yc ~ x), family = categorical(K = 3), data = d))
  expect_identical(variables(fk), variables(fit))
  expect_equal(as.numeric(logLik(fk)), as.numeric(logLik(fit)),
               tolerance = 1e-10)
  # a K that does not match the data keeps the family built from K
  expect_error(frm(bf(yc ~ x), family = categorical(K = 2), data = d))
  d$tr <- 10
  d$y1 <- stats::rbinom(n, 10, 0.3)
  d$y2 <- stats::rbinom(n, 10 - d$y1, 0.5)
  d$y3 <- d$tr - d$y1 - d$y2
  fm <- quiet_fit(frm(bf(cbind(y1, y2, y3) | trials(tr) ~ x),
                      family = multinomial(K = 3), data = d))
  # brms 2.23.0's stancode: muy2, muy3
  expect_identical(names(fixef_by_dpar(fm)), c("muy2", "muy3"))
})
