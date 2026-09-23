# What the recheck of punch round 1 found (dev/reviews, 2026-09-21):
# two silent wrong answers the round-1 fixes introduced, three majors
# and the minors. Every block here was run against the round-1 build
# first and recorded failing in dev/shapes-log/seen-failing-punch2.txt.

p2_data <- local({
  cache <- new.env(parent = emptyenv())
  function() {
    if (!is.null(cache$d)) return(cache$d)
    set.seed(20260921)
    n <- 150
    d <- data.frame(x = rnorm(n), g = factor(rep(1:15, each = 10)))
    u <- rnorm(15, 0, 0.8)
    d$y <- rnorm(n, 1 + 0.5 * d$x + u[d$g], 1)
    d$ord <- factor(cut(0.8 * d$x + u[d$g] + stats::rlogis(n),
                        c(-Inf, -0.5, 0.8, Inf), labels = 1:3),
                    ordered = TRUE)
    e <- matrix(rnorm(2 * n), n, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2))
    d$y1 <- 1 + 0.5 * d$x + e[, 1]
    d$y2 <- -0.3 * d$x + e[, 2]
    d$cnt <- rpois(n, exp(0.3 + 0.4 * d$x))
    cache$d <- d
    d
  }
})

p2_fit <- local({
  cache <- new.env(parent = emptyenv())
  function(key) {
    if (!is.null(cache[[key]])) return(cache[[key]])
    d <- p2_data()
    out <- switch(key,
      mixed = frm(bf(y ~ x + (1 | g)) + gaussian(), data = d),
      ordmix = frm(bf(ord ~ x + (1 | g)) + cumulative(), data = d),
      acatmix = frm(bf(ord ~ x + (1 | g)) + acat(), data = d),
      rescor = frm(bf(mvbind(y1, y2) ~ x) + gaussian() + set_rescor(TRUE),
                   data = d),
      pois = frm(bf(cnt ~ x) + poisson(), data = d),
      ord = frm(bf(ord ~ x) + cumulative(), data = d))
    cache[[key]] <- out
    out
  }
})

# BLOCKER A. Distinct unseen levels draw independent effects.

test_that("each unseen level draws its own effect; one level shares it", {
  fit <- p2_fit("mixed")
  nd <- data.frame(x = 0,
                   g = factor(c("new1", "new2", "new3", "new4", "new4")))
  set.seed(1)
  d <- predict(fit, newdata = nd, allow_new_levels = TRUE,
               propagate_error = FALSE, ndraws = 4000, summary = FALSE)
  r <- stats::cor(d)
  # three DISTINCT unseen levels: brms draws an independent effect for
  # each, so their draws are uncorrelated. They used to share one draw
  # per replicate and correlated at tau^2 / (tau^2 + sigma^2)
  distinct <- r[1:3, 1:3][upper.tri(diag(3))]
  expect_true(all(abs(distinct) < 0.06))
  # two rows in the SAME unseen level load one draw of its effect
  tau2 <- VarCorr(fit)$g$sd[1L, "Estimate"]^2
  s2 <- sigma(fit)^2
  expect_equal(r[4, 5], tau2 / (tau2 + s2), tolerance = 0.06)
})

# BLOCKER D. coef() on a mixed ordinal fit, as brms computes it.

test_that("coef() on a mixed ordinal fit applies brms's threshold rule", {
  # brms:::coef.brmsfit: a random intercept moves each threshold, with
  # the sign the family's specials give it. cumulative and sratio are
  # "thres_minus_eta", so Intercept[k] = threshold_k - r_g; cratio and
  # acat are "eta_minus_thres", so Intercept[k] = r_g - threshold_k.
  # There is no separate intercept column.
  for (key in c("ordmix", "acatmix")) {
    fit <- p2_fit(key)
    cf <- coef(fit)$g
    fe <- fixef(fit)[, "Estimate"]
    r <- ranef(fit)$g[, 1]
    expect_equal(colnames(cf), rownames(fixef(fit)), info = key)
    for (k in 1:2) {
      nm <- paste0("Intercept[", k, "]")
      want <- if (key == "ordmix") fe[[nm]] - r else r - fe[[nm]]
      expect_equal(unname(cf[[nm]]), unname(want), tolerance = 1e-12,
                   info = paste(key, nm))
    }
    # a slope with no group-level term is the same in every group
    expect_equal(length(unique(round(cf$x, 12))), 1L, info = key)
  }
})

# MAJOR B. A row's summary does not depend on its neighbours.

test_that("a row predicted alone and beside an overflowing row agree", {
  fit <- p2_fit("pois")
  b <- fixef_by_dpar(fit)$mu
  # the row where exp(eta) overflows in about half the parameter draws
  x_bad <- (709.5 - b[["(Intercept)"]]) / b[["x"]]
  good <- data.frame(x = 4)
  both <- data.frame(x = c(4, x_bad))
  set.seed(7)
  a <- predict(fit, newdata = good, ndraws = 400)
  set.seed(7)
  w <- NULL
  ab <- withCallingHandlers(predict(fit, newdata = both, ndraws = 400),
    warning = function(cnd) {
      w <<- conditionMessage(cnd)
      invokeRestart("muffleWarning")
    })
  # the overflowing row used to drop the WHOLE replicate, so the good
  # row was summarized over a selected subset of parameter draws
  expect_identical(ab[1L, ], a[1L, ])
  expect_match(w, "row 2")
})

# MAJOR C. Multivariate draws carry the residual correlation.

test_that("predict(summary = FALSE) on a rescor fit draws jointly", {
  fit <- p2_fit("rescor")
  rc <- rescor_matrix(fit)[1L, 2L]
  set.seed(3)
  d <- predict(fit, summary = FALSE, ndraws = 3000,
               propagate_error = FALSE)
  within <- vapply(seq_len(dim(d)[2L]), function(i) {
    stats::cor(d[, i, 1L], d[, i, 2L])
  }, 0)
  # it drew each response on its own and correlated at 0.009
  expect_equal(mean(within), rc, tolerance = 0.03)
})

# Minors.

test_that("ranef() and coef() both say Intercept, as brms does", {
  fit <- p2_fit("mixed")
  expect_equal(colnames(ranef(fit)$g), "Intercept")
  expect_true("Intercept" %in% colnames(coef(fit)$g))
})

test_that("insight::get_parameters() reports the thresholds fixef() does", {
  skip_if_not_installed("insight")
  fit <- p2_fit("ord")
  p <- insight::get_parameters(fit)
  fe <- fixef(fit)
  expect_true(all(c("Intercept[1]", "Intercept[2]") %in% p$Parameter))
  expect_equal(p$Estimate[match(rownames(fe), p$Parameter)],
               unname(fe[, "Estimate"]))
  V <- insight::get_varcov(fit)
  expect_equal(rownames(V), as.character(p$Parameter))
})

test_that("predict(summary = FALSE) names no draw and no row, as brms", {
  # brms's predict(summary = FALSE) and posterior_predict() on the same
  # kind of fit return dimnames list(NULL, NULL), and list(NULL, NULL,
  # c("y1", "y2")) on a multivariate one (dev/shapes-p2-brmsref.R). The
  # draws kept the data's row names on the observation margin while the
  # summary dropped them.
  fit <- p2_fit("pois")
  set.seed(4)
  d <- predict(fit, summary = FALSE, ndraws = 5)
  expect_true(all(vapply(dimnames(d) %||% list(NULL), is.null, TRUE)))
  set.seed(4)
  m <- predict(p2_fit("rescor"), summary = FALSE, ndraws = 5)
  expect_null(dimnames(m)[[1L]])
  expect_null(dimnames(m)[[2L]])
  expect_equal(dimnames(m)[[3L]], c("y1", "y2"))
})
