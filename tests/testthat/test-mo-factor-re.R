# mo() crossed with a factor, and mo() in group-level terms
# (dev/moparity-findings.md, item 1). Both were refused on 0.69.0
# (rellib-r7, dev/moparity-refusals.R): "mo() interactions support
# numeric multipliers only: f is factor" and "mo() is only supported as
# a standalone additive term or a two-way ':'/'*' interaction".
#
# brms 2.23.0 writes wrong Stan code for a factor crossed with mo()
# (brms#1828, acknowledged): every column is multiplied by the first
# covariate Csp_1, and a factor of three or more levels stops the model
# build (dev/moparity-brms-code.R). The model fitted here is the one
# its Stan data describe and brms 2.22.0 fits: one simplex and one
# coefficient per coded column. The log densities below are checked
# against brms's program for the same model written with numeric
# indicator columns, where brms 2.23.0 is right.

mo_fac_data <- function(seed = 1, n = 300, n_g = 20) {
  set.seed(seed)
  dd <- data.frame(
    x = sample(0:3, n, TRUE),
    f = factor(sample(c("a", "b", "c"), n, TRUE)),
    g = factor(sample(paste0("g", seq_len(n_g)), n, TRUE)),
    z = stats::rnorm(n), w = stats::rnorm(n)
  )
  eff <- c(0, 1, 1.5, 1.7)
  bg <- stats::rnorm(n_g, 0, 0.5)
  ag <- stats::rnorm(n_g, 0, 0.4)
  dd$y <- 1 + ag[dd$g] +
    (1 + (dd$f == "b") - (dd$f == "c") + bg[dd$g]) * eff[dd$x + 1] +
    stats::rnorm(n, 0, 0.5)
  # the same factor as numeric indicators, the spelling brms 2.23.0
  # writes correct code for
  dd$fa <- as.numeric(dd$f == "a")
  dd$fb <- as.numeric(dd$f == "b")
  dd$fc <- as.numeric(dd$f == "c")
  dd
}

test_that("mo(x) * f codes one simplex per contrast column", {
  dd <- mo_fac_data()
  fit <- frm(bf(y ~ mo(x) * f), family = gaussian(), data = dd)
  expect_identical(rownames(fixef(fit)),
                   c("Intercept", "fb", "fc", "mox", "mox:fb", "mox:fc"))
  expect_identical(grep("^zeta", names(fit$estimates), value = TRUE),
                   c("zeta1", "zeta2", "zeta3"))
  expect_true(all(c("bsp_mox", "bsp_mox:fb", "bsp_mox:fc") %in%
                    variables(fit)))
  # the same model written with numeric indicators
  ind <- frm(bf(y ~ fb + fc + mo(x) + mo(x):fb + mo(x):fc),
             family = gaussian(), data = dd)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ind)),
               tolerance = 1e-8)
  expect_equal(unname(fixef(fit)[, 1]), unname(fixef(ind)[, 1]),
               tolerance = 1e-5)
})

test_that("mo(x):f without the main term codes every level", {
  dd <- mo_fac_data()
  fit <- frm(bf(y ~ mo(x):f), family = gaussian(), data = dd)
  expect_identical(rownames(fixef(fit)),
                   c("Intercept", "mox:fa", "mox:fb", "mox:fc"))
  ind <- frm(bf(y ~ mo(x):fa + mo(x):fb + mo(x):fc), family = gaussian(),
             data = dd)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ind)),
               tolerance = 1e-8)
  # brms's name order follows terms(): f is read first here
  fit2 <- frm(bf(y ~ f + mo(x):f), family = gaussian(), data = dd)
  expect_identical(rownames(fixef(fit2))[4:6],
                   c("fa:mox", "fb:mox", "fc:mox"))
  fit3 <- frm(bf(y ~ z + mo(x):z), family = gaussian(), data = dd)
  expect_identical(rownames(fixef(fit3))[3], "z:mox")
})

test_that("a three-way mo() interaction is brms's sp columns", {
  dd <- mo_fac_data()
  fit <- frm(bf(y ~ mo(x) * z * w), family = gaussian(), data = dd)
  expect_identical(rownames(fixef(fit)),
                   c("Intercept", "z", "w", "z:w", "mox", "mox:z", "mox:w",
                     "mox:z:w"))
  # the product column is z * w, as model.matrix() codes it
  lp <- fit$frame$linpreds[[1L]]
  m3 <- Filter(function(m) identical(m$label, "mox:z:w"), lp$mo)[[1L]]
  expect_identical(m3$mult, dd$z * dd$w)
})

test_that("a factor mo() interaction predicts new rows as in sample", {
  dd <- mo_fac_data()
  fit <- frm(bf(y ~ mo(x) * f), family = gaussian(), data = dd)
  nd <- dd[c(5, 1, 9, 2), ]
  expect_equal(fitted(fit, newdata = nd)[, 1], fitted(fit)[c(5, 1, 9, 2), 1],
               tolerance = 1e-12)
  # one row alone keeps the fitted levels and contrasts
  expect_equal(fitted(fit, newdata = dd[2, ])[, 1], fitted(fit)[2, 1],
               tolerance = 1e-12)
  bad <- dd[1, ]
  bad$f <- factor("d")
  expect_error(fitted(fit, newdata = bad), "d")
})

test_that("a group-level mo() term varies the population coefficient", {
  dd <- mo_fac_data()
  fit <- allow_warnings(
    frm(bf(y ~ mo(x) + (mo(x) | g)), family = gaussian(), data = dd),
    "Boundary")
  vc <- VarCorr(fit)$g$sd
  expect_identical(rownames(vc), c("Intercept", "mox"))
  expect_true(all(c("sd_g__mox", "cor_g__Intercept__mox") %in%
                    variables(fit)))
  # the in-sample linear predictor is (bsp + r_g) * mo(x) + r0_g
  w <- frmtmb:::mo_simplex(fit$estimates$zeta1)
  mo_x <- 3 * c(0, cumsum(w))[dd$x + 1]
  re <- ranef(fit)$g[as.character(dd$g), , drop = FALSE]
  eta <- fixef(fit)["Intercept", 1] + re[, "Intercept"] +
    (fixef(fit)["mox", 1] + re[, "mox"]) * mo_x
  expect_equal(unname(fitted(fit)[, 1]), unname(eta), tolerance = 1e-10)
  nd <- dd[c(7, 3, 11), ]
  expect_equal(fitted(fit, newdata = nd)[, 1], fitted(fit)[c(7, 3, 11), 1],
               tolerance = 1e-12)
  # population level only
  expect_equal(unname(fitted(fit, newdata = nd, re_formula = NA)[, 1]),
               unname(fixef(fit)["Intercept", 1] +
                        fixef(fit)["mox", 1] * mo_x[c(7, 3, 11)]),
               tolerance = 1e-10)
})

test_that("a group-level mo() term needs its population-level column", {
  dd <- mo_fac_data()
  # brms 2.23.0 drops (mo(x) | g) without a word when the formula has
  # no population-level mo(x): its Stan code reads no simplex for it
  expect_error(frm(bf(y ~ (mo(x) | g)), family = gaussian(), data = dd),
               "require corresponding population-level terms")
  expect_error(frm(bf(y ~ mo(x) + (mo(x):z | g)), family = gaussian(),
                   data = dd),
               "require corresponding population-level terms")
  # the guard absent: with the column present the same term fits
  fit <- allow_warnings(frm(bf(y ~ mo(x) * z + (0 + mo(x):z | g)),
                            family = gaussian(), data = dd),
                        c("Boundary", "singular"))
  expect_s3_class(fit, "frmtmb_fit")
  expect_error(frm(bf(y ~ mo(x) + (1 | mo(x))), family = gaussian(),
                   data = dd), "grouping factor")
  expect_error(frm(bf(y ~ mo(x) + (mo(x) | g)), family = gaussian(),
                   data = dd, quadrature = TRUE), "quadrature")
})

# ---------------------------------------------------------------------
# Against brms's compiled programs (gated: Stan compiles).
# ---------------------------------------------------------------------

test_that("mo(x) * f is brms's model written with indicators", {
  skip_unless_brms_fit()
  dd <- mo_fac_data()
  fit <- frm(bf(y ~ mo(x) * f), family = gaussian(), data = dd)
  # three simplexes, each with brms's flat Dirichlet constant lgamma(3)
  brms_lp_check(brms::bf(y ~ fb + fc + mo(x) + mo(x):fb + mo(x):fc),
                gaussian(), dd, fit, const = 3 * lgamma(3))
  fit2 <- frm(bf(y ~ mo(x):f), family = gaussian(), data = dd)
  brms_lp_check(brms::bf(y ~ mo(x):fa + mo(x):fb + mo(x):fc), gaussian(),
                dd, fit2, const = 3 * lgamma(3))
})

test_that("mo() names match brms's for a factor and a three-way term", {
  skip_unless_brms()
  dd <- mo_fac_data()
  for (fo in list(y ~ mo(x) * f, y ~ mo(x):f, y ~ f + mo(x):f,
                  y ~ z + mo(x):z, y ~ mo(x) * z * w)) {
    fit <- frm(bf(fo), family = gaussian(), data = dd)
    lp <- fit$frame$linpreds[[1L]]
    ours <- vapply(lp$mo, `[[`, "", "label")
    # brms 2.23.0 builds its table for these from the sp model matrix,
    # which is right; only its Stan code is not (brms#1828). A factor of
    # three levels crossed with the main term stops that table too, so
    # the two-level case stands in for it
    dfo <- if (identical(fo, y ~ mo(x) * f)) y ~ mo(x) * f2 else fo
    dd2 <- dd
    dd2$f2 <- factor(ifelse(dd$f == "a", "a", "b"))
    dp <- suppressMessages(brms::default_prior(dfo, data = dd2))
    simo <- sort(dp$coef[dp$class == "simo"])
    if (identical(fo, y ~ mo(x) * f)) {
      fit2 <- frm(bf(dfo), family = gaussian(), data = dd2)
      ours <- vapply(fit2$frame$linpreds[[1L]]$mo, `[[`, "", "label")
    }
    expect_identical(sort(paste0(ours, "1")), simo,
                     info = deparse1(fo))
  }
})

test_that("a group-level mo() term is brms's (bsp + r) * mo(simo, x)", {
  skip_unless_brms_fit()
  dd <- mo_fac_data(seed = 4, n = 400, n_g = 25)
  fit <- frm(bf(y ~ mo(x) + (mo(x) | g)), family = gaussian(), data = dd)
  brms_lp_check(brms::bf(y ~ mo(x) + (mo(x) | g)), gaussian(), dd, fit,
                joint = TRUE, const = lgamma(3))
  fit2 <- allow_warnings(
    frm(bf(y ~ mo(x) * f + (1 + mo(x) || g)), family = gaussian(),
        data = dd), "Boundary")
  brms_lp_check(brms::bf(y ~ fb + fc + mo(x) + mo(x):fb + mo(x):fc +
                           (1 + mo(x) || g)),
                gaussian(), dd, fit2, joint = TRUE, const = 3 * lgamma(3))
})

# Punch round 1, m7: a `||` term with mo() crossed with a factor.
# Seen to fail on the lane build before the fix: lme4's piecewise
# expansion coded `0 + mo(x):f | g` with every level of f, and the term
# was refused for "mox:fa", a population-level term it did not need.
# brms reads `||` as one block with no correlations, coded whole
# (Intercept, fb, fc, mox, mox:fb, mox:fc; dev/moparity-rev-dblbar.R),
# and so does frmtmb when the term has mo().
test_that("(mo(x) * f || g) is one uncorrelated block, coded whole", {
  dd <- mo_fac_data()
  fit <- allow_warnings(
    frm(bf(y ~ mo(x) * f + (mo(x) * f || g)), family = gaussian(),
        data = dd),
    c("Boundary", "singular"))
  v <- variables(fit)
  expect_setequal(grep("^sd_g__", v, value = TRUE),
                  paste0("sd_g__", c("Intercept", "fb", "fc", "mox",
                                     "mox:fb", "mox:fc")))
  expect_identical(vapply(fit$frame$re_blocks, `[[`, "", "covstruct"),
                   "diag")
  # the factor written as indicators is the same model: the coding is
  # brms's, which the indicator spelling states column by column
  ind <- allow_warnings(
    frm(bf(y ~ fb + fc + mo(x) + mo(x):fb + mo(x):fc +
             (fb + fc + mo(x) + mo(x):fb + mo(x):fc || g)),
        family = gaussian(), data = dd),
    c("Boundary", "singular"))
  ll <- as.numeric(logLik(fit))
  expect_lt(abs(ll - as.numeric(logLik(ind))), 1e-8 * abs(ll))
})

# Punch round 1, m8. Seen to fail on the lane build before the fix: the
# refusal said every block must be a dim-1 us/diag term, which
# (0 + mo(x) | g) is.
test_that("quadrature names group-level mo() as the reason it refuses", {
  dd <- mo_fac_data()
  expect_error(frm(bf(y ~ mo(x) + (0 + mo(x) | g)), family = gaussian(),
                   data = dd, quadrature = TRUE),
               "does not support group-level mo")
  # the guard absent: a population-level mo() with a scalar intercept
  fit <- allow_warnings(frm(bf(y ~ mo(x) + (1 | g)), family = gaussian(),
                            data = dd, quadrature = TRUE),
                        c("Boundary", "singular"))
  expect_s3_class(fit, "frmtmb_fit")
})
