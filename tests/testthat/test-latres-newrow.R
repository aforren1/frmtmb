# New rows of a latent residual, beside what else the model holds. Each block here
# read a wrong number in silence on the first build of 0.79.0, where seven lanes
# met (dev/reviews/2026-11-04-release.md, findings 1, 2, 3, 8, 9, 13). The references
# hold no frmtmb reader: a fit with no latent, the draw's own closed form.

lrn_data <- function(seed, sd_e = 0.06) {
  set.seed(seed)
  ng <- 8
  nt <- 25
  d <- expand.grid(t = 1:nt, s = factor(1:ng))
  d$x <- stats::rnorm(nrow(d))
  re <- stats::rnorm(ng, 0, 0.8)
  e <- stats::ave(stats::rnorm(nrow(d)), d$s, FUN = function(z) {
    stats::filter(z, 0.5, "recursive")
  })
  d$e <- e
  d$y <- stats::rpois(nrow(d), exp(3 + 0.3 * d$x + re[d$s] + sd_e * e))
  d
}

test_that("a seen level's Est.Error beside a latent holds no more than the two parts", {
  # the differences moved beta and theta at a fixed b, so the intercept's whole
  # variance reached the row: 8.35, where the reading of the two parts is 1.77
  for (seed in c(9101, 7311, 7312)) {
    d <- lrn_data(seed)
    nd <- data.frame(t = 26, s = factor(c(1, 2)), x = c(0, 1))
    f <- suppressWarnings(frm(y ~ x + (1 | s) + ar(time = t, gr = s, cov = TRUE),
                              data = d, family = poisson()))
    f0 <- suppressWarnings(frm(y ~ x + (1 | s), data = d, family = poisson()))
    sl <- latres_row_sd(f, "y", NULL)
    a <- fitted(f, newdata = nd)
    b <- fitted(f0, newdata = nd)
    ref <- sqrt(b[, "Est.Error"]^2 + (a[, "Estimate"] * sqrt(exp(sl^2) - 1))^2)
    expect_lt(max(abs(a[, "Est.Error"] / ref - 1)), 0.3, label = paste("seed", seed))
    expect_lt(max(a[, "Q97.5"] / a[, "Estimate"]), 1.35, label = paste("seed", seed))
  }
})

test_that("a field's far row beside a latent keeps its kriging", {
  # at a position far from the data the field's prior variance is the row's, and the
  # differences at a fixed z read none of it
  for (seed in c(402, 7313)) {
    set.seed(seed)
    n <- 100
    d <- data.frame(t = 1:n, s = factor(1), x = sort(stats::runif(n, 0, 10)))
    d$y <- stats::rpois(n, exp(3 + sin(d$x) + 0.2 * as.numeric(
      stats::filter(stats::rnorm(n), 0.5, "recursive"))))
    f <- suppressWarnings(frm(y ~ gp(x) + ar(time = t, gr = s, cov = TRUE),
                              data = d, family = poisson()))
    f0 <- suppressWarnings(frm(y ~ gp(x), data = d, family = poisson()))
    nd <- data.frame(t = n + 1L, s = factor(1), x = 30)
    a <- fitted(f, newdata = nd)
    b <- fitted(f0, newdata = nd)
    # no narrower than the no-latent read, to 10 percent: the kriging is there
    expect_gt(a[, "Est.Error"], 0.9 * b[, "Est.Error"], label = paste("seed", seed))
  }
})

test_that("predictions() reads each response's own mean, in either order", {
  # with no resp the frame holds every response, and only the first was read: a
  # me() or a latent on another response read the plug-in, 14 percent low
  skip_if_not_installed("marginaleffects")
  for (seed in c(9102, 7321)) {
    set.seed(seed)
    n <- 200
    d <- data.frame(x = stats::rnorm(n), zo = stats::rnorm(n), sz = 0.6)
    d$yg <- stats::rnorm(n, 1 + d$x)
    d$yq <- stats::rpois(n, exp(0.5 + 0.5 * d$x + 0.6 * d$zo))
    nd <- data.frame(x = c(-1, 0, 1), zo = 0.5, sz = 0.6)
    both <- list(
      list(mvbf(bf(yg ~ x), bf(yq ~ x + me(zo, sz))), list(gaussian(), poisson())),
      list(mvbf(bf(yq ~ x + me(zo, sz)), bf(yg ~ x)), list(poisson(), gaussian())))
    for (one in both) {
      fit <- suppressWarnings(frm(one[[1L]], family = one[[2L]], data = d))
      ft <- fitted(fit, newdata = nd)
      p <- suppressWarnings(marginaleffects::predictions(fit, newdata = nd, vcov = FALSE))
      for (r in c("yg", "yq")) {
        got <- p$estimate[as.character(p$group) == r]
        expect_lt(max(abs(got / ft[, "Estimate", r] - 1)), 1e-8,
                  label = paste("seed", seed, r))
      }
    }
  }
})

test_that("predictions() on a latent on the second response reads its mean", {
  skip_if_not_installed("marginaleffects")
  set.seed(7322)
  ng <- 8
  nt <- 25
  d <- expand.grid(t = 1:nt, s = factor(1:ng))
  d$x <- stats::rnorm(nrow(d))
  e <- stats::ave(stats::rnorm(nrow(d)), d$s, FUN = function(z) {
    stats::filter(z, 0.5, "recursive")
  })
  d$yg <- stats::rnorm(nrow(d), 1 + d$x)
  d$yp <- stats::rpois(nrow(d), exp(0.5 + 0.4 * d$x + 0.5 * e))
  fit <- suppressWarnings(frm(mvbf(bf(yg ~ x), bf(yp ~ x + ar(time = t, gr = s, cov = TRUE))),
                              family = list(gaussian(), poisson()), data = d))
  nd <- data.frame(t = 26, s = factor(1), x = c(-1, 0, 1))
  ft <- fitted(fit, newdata = nd)
  p <- suppressWarnings(marginaleffects::predictions(fit, newdata = nd, vcov = FALSE))
  got <- p$estimate[as.character(p$group) == "yp"]
  expect_lt(max(abs(got / ft[, "Estimate", "yp"] - 1)), 1e-8)
})

test_that("a me() on shape beside a latent on mu keeps the latent", {
  # the mean was low by exp(s^2 / 2), Est.Error up to 7.6 times narrow, and the
  # interval had zero width
  for (seed in c(9103, 7331)) {
    set.seed(seed)
    ng <- 8
    nt <- 25
    d <- expand.grid(t = 1:nt, s = factor(1:ng))
    d$x <- stats::rnorm(nrow(d))
    d$zo <- stats::rnorm(nrow(d))
    d$sz <- 0.3
    e <- stats::ave(stats::rnorm(nrow(d)), d$s, FUN = function(z) {
      stats::filter(z, 0.5, "recursive")
    })
    d$y <- stats::rgamma(nrow(d), shape = 5, rate = 5 / exp(0.5 + 0.4 * d$x + 0.4 * e))
    f <- suppressWarnings(frm(bf(y ~ x + ar(time = t, gr = s, cov = TRUE), shape ~ me(zo, sz)),
                              data = d, family = Gamma(link = "log")))
    nd <- data.frame(t = 26, s = factor(c(1, 2)), x = c(0, 1), zo = 0, sz = 0.3)
    sl <- latres_row_sd(f, "y", NULL)
    eta <- as.numeric(frm_linpred(f, newdata = nd, dpar = "mu", type = "link"))
    a <- fitted(f, newdata = nd)
    expect_lt(max(abs(a[, "Estimate"] / exp(eta + sl^2 / 2) - 1)), 1e-6,
              label = paste("seed", seed))
    # the draw's own spread
    expect_gt(min(a[, "Est.Error"] / (exp(eta + sl^2 / 2) * sqrt(exp(sl^2) - 1))), 0.9,
              label = paste("seed", seed))
    expect_gt(min(a[, "Q97.5"] - a[, "Q2.5"]), 0, label = paste("seed", seed))
  }
})

test_that("a category's Est.Error beside a gp() and a latent keeps the kriging", {
  # the category reader replaced eta_shift, which held fit_fd_se()'s chain step, so
  # d out / d eta read 0
  for (seed in c(1501, 7332)) {
    set.seed(seed)
    n <- 120
    d <- data.frame(t = 1:n, s = factor(1), pos = sort(stats::runif(n, 0, 10)))
    e <- as.numeric(stats::filter(stats::rnorm(n), 0.5, "recursive"))
    d$o <- ordered(cut(sin(d$pos) + 1.5 * e + stats::rlogis(n), c(-Inf, -1, 0, 1, Inf),
                       labels = 1:4))
    f <- suppressWarnings(frm(o ~ gp(pos) + ar(time = t, gr = s, cov = TRUE),
                              data = d, family = acat()))
    f0 <- suppressWarnings(frm(o ~ gp(pos), data = d, family = acat()))
    nd <- data.frame(t = n + 1L, s = factor(1), pos = 30)
    a <- suppressWarnings(fitted(f, newdata = nd))
    b <- suppressWarnings(fitted(f0, newdata = nd))
    # category 1 at the far row
    expect_gt(a[1L, "Est.Error", 1L], 0.9 * b[1L, "Est.Error", 1L],
              label = paste("seed", seed))
  }
})

test_that("the linear scale holds the latent's draw, as a me() row does", {
  d <- lrn_data(7333, 0.3)
  f <- suppressWarnings(frm(y ~ x + (1 | s) + ar(time = t, gr = s, cov = TRUE),
                            data = d, family = poisson()))
  nd <- data.frame(t = 26, s = factor(c(1, 2)), x = c(0, 1))
  sl <- latres_row_sd(f, "y", NULL)
  se0 <- as.numeric(suppressWarnings(frm_linpred(f, newdata = nd, type = "link",
                                                  se.fit = TRUE))$se.fit)
  a <- fitted(f, newdata = nd, scale = "linear")
  expect_equal(unname(a[, "Est.Error"]), sqrt(se0^2 + sl^2), tolerance = 1e-8)
})

test_that("a standard error refused above the dense limit says so", {
  d <- lrn_data(7334, 0.3)
  f <- suppressWarnings(frm(y ~ x + (1 | s) + ar(time = t, gr = s, cov = TRUE),
                            data = d, family = poisson()))
  nd <- data.frame(t = 26, s = factor(c(1, 2)), x = c(0, 1))
  withr::local_options(frmtmb.latres_dense_rows = 10L, frmtmb.latres_se_rows = 5L)
  said <- 0L
  count <- function(m) {
    said <<- said + 1L
    invokeRestart("muffleMessage")
  }
  v <- withCallingHandlers(fitted(f, newdata = nd, scale = "linear"),
                           frmtmb_latres_se_skipped = count)
  expect_gte(said, 1L)
  expect_true(all(is.finite(v[, "Estimate"])))
  said <- 0L
  v <- withCallingHandlers(fitted(f), frmtmb_latres_se_skipped = count)
  expect_gte(said, 1L)
})

test_that("new rows above the row limit beside another block read the joint covariance", {
  d <- lrn_data(7335, 0.3)
  fm <- y ~ x + (1 | s) + ar(time = t, gr = s, cov = TRUE)
  f <- suppressWarnings(frm(fm, data = d, family = poisson()))
  nd <- data.frame(t = 26, s = factor(c(1, 2)), x = c(0, 1))
  a <- fitted(f, newdata = nd)
  withr::with_options(list(frmtmb.latres_se_rows = 5L), {
    expect_equal(fitted(f, newdata = nd)[, "Est.Error"], a[, "Est.Error"], tolerance = 1e-8)
  })
  # past the dense limit it reads NA and says so; a new fit, since the matrix is cached
  f2 <- suppressWarnings(frm(fm, data = d, family = poisson()))
  withr::local_options(frmtmb.latres_se_rows = 5L, frmtmb.latres_dense_rows = 10L)
  said <- 0L
  v <- withCallingHandlers(fitted(f2, newdata = nd),
                           frmtmb_latres_se_skipped = function(m) {
                             said <<- said + 1L
                             invokeRestart("muffleMessage")
                           })
  expect_gte(said, 1L)
  expect_true(all(is.na(v[, "Est.Error"])))
  expect_true(all(is.finite(v[, "Estimate"])))
})

test_that("above the row limit, Est.Error beside a field or a smooth holds its two parts", {
  # real rows, and a latent small enough that the term's own part shows: the mean's
  # closed form, and the no-latent fit's Est.Error with the draw in quadrature
  set.seed(7552)
  n <- 3001L
  nt <- 50L
  ng <- ceiling(n / nt)
  d <- expand.grid(t = seq_len(nt), s = factor(seq_len(ng)))[seq_len(n), ]
  d$x <- stats::runif(n, -2, 2)
  e <- stats::ave(stats::rnorm(n), d$s, FUN = function(z) stats::filter(z, 0.5, "recursive"))
  d$y <- stats::rpois(n, exp(1 + sin(d$x) + 0.04 * e))
  nd <- data.frame(t = nt + 1L, s = factor(c(1, 2), levels = levels(d$s)), x = c(0, 1.5))
  for (term in c("gp(x, k = 10)", "s(x)")) {
    f <- suppressWarnings(frm(stats::as.formula(
      paste0("y ~ ", term, " + ar(time = t, gr = s, cov = TRUE)")), data = d, family = poisson()))
    f0 <- suppressWarnings(frm(stats::as.formula(paste0("y ~ ", term)), data = d, family = poisson()))
    sl <- latres_row_sd(f, "y", NULL)
    a <- fitted(f, newdata = nd)
    b <- fitted(f0, newdata = nd)
    eta <- as.numeric(frm_linpred(f, newdata = nd, type = "link", dpar = "mu"))
    expect_lt(max(abs(a[, "Estimate"] / exp(eta + sl^2 / 2) - 1)), 1e-6, label = term)
    two <- sqrt(b[, "Est.Error"]^2 + (a[, "Estimate"] * sqrt(exp(sl^2) - 1))^2)
    expect_lt(max(abs(a[, "Est.Error"] / two - 1)), 0.05, label = term)
  }
})
