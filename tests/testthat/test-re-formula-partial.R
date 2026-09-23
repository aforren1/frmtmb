# A one-sided re_formula keeps the group-level terms it names and drops
# the rest, as brms's update_re_terms() does. Until this lane it was a
# two-way switch: a formula naming SOME terms returned the prediction
# with ALL of them (dev/adefects-findings.md section 11, item 7), and a
# formula naming a grouping factor the fit does not have did the same
# (dev/test-backlog.md). dev/reunc-brms.R measures the same selections
# against brms 2.23.0 on one data set.

partial_fixture <- function() {
  set.seed(3)
  G <- 15
  m <- 8
  d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                  h = factor(rep(1:5, times = G * m / 5)),
                  x = stats::rnorm(G * m))
  ug <- stats::rnorm(G, 0, 0.8)
  sg <- stats::rnorm(G, 0, 0.4)
  uh <- stats::rnorm(5, 0, 0.6)
  d$y <- 1 + 0.5 * d$x + ug[d$g] + sg[d$g] * d$x + uh[d$h] +
    stats::rnorm(G * m)
  d
}

partial_fit <- function() {
  frm(bf(y ~ x + (1 + x | g) + (1 | h)), data = partial_fixture())
}

# the linear predictor built by hand from fixef() and ranef(): which
# columns of which term enter is the whole question
partial_hand <- function(fit, d, ig, sl, ih) {
  fe <- fixef(fit)[, "Estimate"]
  re <- ranef(fit)
  fe[[1]] + fe[[2]] * d$x +
    ig * re$g[as.character(d$g), "Intercept"] +
    sl * re$g[as.character(d$g), "x"] * d$x +
    ih * re$h[as.character(d$h), "Intercept"]
}

test_that("a partial re_formula keeps exactly the terms it names", {
  fit <- partial_fit()
  d <- partial_fixture()
  cases <- list(
    list(f = ~(1 | g), k = c(1, 0, 0)),
    list(f = ~(0 + x | g), k = c(0, 1, 0)),
    list(f = ~(1 + x | g), k = c(1, 1, 0)),
    list(f = ~(1 | h), k = c(0, 0, 1)),
    list(f = ~(1 | g) + (1 | h), k = c(1, 0, 1))
  )
  for (cs in cases) {
    want <- partial_hand(fit, d, cs$k[1], cs$k[2], cs$k[3])
    expect_equal(unname(frm_linpred(fit, re_formula = cs$f)), unname(want),
                 info = deparse1(cs$f))
    expect_equal(unname(fitted(fit, re_formula = cs$f)[, "Estimate"]),
                 unname(want), info = deparse1(cs$f))
  }
  # every column of every term is the full prediction, bit for bit
  expect_identical(frm_linpred(fit, re_formula = ~(1 + x | g) + (1 | h)),
                   frm_linpred(fit))
  expect_identical(fitted(fit, re_formula = ~(x | g) + (1 | h)),
                   fitted(fit))
})

test_that("a formula with no group-level term is the population level", {
  fit <- partial_fit()
  # brms reads ~1 and ~0 as "no group-level effects"; ~1 used to keep
  # every term
  expect_identical(frm_linpred(fit, re_formula = ~1),
                   frm_linpred(fit, re_formula = NA))
  expect_identical(fitted(fit, re_formula = ~0), fitted(fit, re_formula = NA))
})

test_that("a dropped term's grouping column is not needed in newdata", {
  fit <- partial_fit()
  d <- partial_fixture()
  rows <- c(1, 9, 17)
  nd <- d[rows, c("x", "g")]
  want <- partial_hand(fit, d, 1, 0, 0)[rows]
  expect_equal(unname(frm_linpred(fit, newdata = nd, re_formula = ~(1 | g))),
               unname(want))
  # and a level of the dropped factor the fit never saw is not refused
  nd$h <- factor("never")
  expect_equal(unname(fitted(fit, newdata = nd,
                             re_formula = ~(1 | g))[, "Estimate"]),
               unname(want))
  expect_error(frm_linpred(fit, newdata = nd), class = "frmtmb_error")
})

test_that("only the kept terms contribute their uncertainty", {
  fit <- partial_fit()
  d <- partial_fixture()
  # the delta method over the kept columns, by hand: the fixed effects
  # and the h block, from the joint covariance
  jc <- frm_joint_cov(fit)
  rn <- jc$names
  fit_h <- frm_linpred(fit, re_formula = ~(1 | h), se.fit = TRUE)
  bk <- fit$frame[["re_blocks"]][[2]]
  X <- cbind(1, d$x)
  Zh <- stats::model.matrix(~ 0 + h, d)
  A <- cbind(X, Zh)
  pos <- c(which(rn == "beta"), which(rn == "b")[bk[["b_idx"]]])
  V <- as.matrix(jc$V[pos, pos])
  expect_equal(unname(fit_h$se.fit), unname(sqrt(rowSums((A %*% V) * A))))
  # and a smaller set is a different number from the full set
  full <- frm_linpred(fit, se.fit = TRUE)$se.fit
  expect_false(isTRUE(all.equal(unname(full), unname(fit_h$se.fit))))
})

test_that("predict() honors a partial re_formula", {
  fit <- partial_fit()
  d <- partial_fixture()
  nd <- d[c(1, 9, 17), c("x", "g")]
  set.seed(4)
  P <- predict(fit, newdata = nd, re_formula = ~(1 | g), summary = FALSE,
               ndraws = 2000)
  want <- fitted(fit, newdata = nd, re_formula = ~(1 | g))[, "Estimate"]
  z <- (colMeans(P) - want) / (apply(P, 2, stats::sd) / sqrt(nrow(P)))
  expect_true(all(abs(z) < 4), info = paste(round(z, 2), collapse = " "))
})

test_that("a term the fit does not have is refused, naming it", {
  fit <- partial_fit()
  for (f in list(~(1 | nosuch), ~(1 + z | g), ~(1 | g) + (1 | nosuch))) {
    expect_error(frm_linpred(fit, re_formula = f), class = "frmtmb_error",
                 regexp = "matches no group-level term")
    expect_error(fitted(fit, re_formula = f), class = "frmtmb_error",
                 regexp = "matches no group-level term")
    expect_error(predict(fit, re_formula = f, ndraws = 5),
                 class = "frmtmb_error", regexp = "matches no group-level")
  }
  e <- tryCatch(predict(fit, re_formula = ~(1 | nosuch), ndraws = 5),
                error = function(e) conditionMessage(e))
  expect_match(e, "(1 | nosuch)", fixed = TRUE)
  expect_match(e, "(1 + x | g)", fixed = TRUE)
})

test_that("brms's spellings of a term resolve", {
  fit <- partial_fit()
  keep_g <- frm_linpred(fit, re_formula = ~(1 | g))
  # an id, a gr() wrapper and a population term beside the bar change
  # nothing about which columns are kept
  expect_identical(frm_linpred(fit, re_formula = ~(1 | p | g)), keep_g)
  expect_identical(frm_linpred(fit, re_formula = ~(1 | gr(g))), keep_g)
  expect_identical(frm_linpred(fit, re_formula = ~ x + (1 | g)), keep_g)
})

test_that("a term shared by two distributional parameters is kept in both", {
  set.seed(8)
  d <- data.frame(g = factor(rep(1:12, each = 10)),
                  h = factor(rep(1:6, 20)), x = stats::rnorm(120))
  d$y <- 1 + d$x + stats::rnorm(12)[d$g] + stats::rnorm(6, 0, 0.5)[d$h] +
    stats::rnorm(120, 0, exp(stats::rnorm(12, 0, 0.3)[d$g]))
  fit <- frm(bf(y ~ x + (1 | g) + (1 | h), sigma ~ (1 | g)), data = d)
  s_null <- frm_linpred(fit, dpar = "sigma")
  s_g <- frm_linpred(fit, dpar = "sigma", re_formula = ~(1 | g))
  s_h <- frm_linpred(fit, dpar = "sigma", re_formula = ~(1 | h))
  s_na <- frm_linpred(fit, dpar = "sigma", re_formula = NA)
  expect_identical(s_g, s_null)
  expect_identical(s_h, s_na)
  mu_g <- frm_linpred(fit, re_formula = ~(1 | g))
  expect_false(isTRUE(all.equal(mu_g, frm_linpred(fit))))
})

test_that("a partial formula is refused where a term cannot be named", {
  set.seed(9)
  d <- data.frame(g = factor(rep(1:6, each = 20)), h = factor(rep(1:4, 30)),
                  k = factor(rep(1:5, each = 4, length.out = 120)),
                  x = stats::runif(120))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(6, 0, 0.3)[d$g] +
    stats::rnorm(4, 0, 0.3)[d$h] + stats::rnorm(5, 0, 0.3)[d$k] +
    stats::rnorm(120, 0, 0.3)
  fit <- suppressWarnings(
    frm(bf(y ~ s(x, g, bs = "fs", k = 4) + (1 | h) + (1 | k)), data = d))
  # the factor smooth is group-level content re_formula = NA drops, and
  # a formula has no way to name it
  expect_error(frm_linpred(fit, re_formula = ~(1 | h)),
               class = "frmtmb_error", regexp = "cannot name")
  # naming every bar term keeps the smooth, as NULL does
  expect_identical(frm_linpred(fit, re_formula = ~(1 | h) + (1 | k)),
                   frm_linpred(fit))
})
