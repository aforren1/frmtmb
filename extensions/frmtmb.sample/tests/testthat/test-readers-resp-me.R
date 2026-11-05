# The readers of the fitted predictor on draws: me() on new data, cs() on the
# linear scale, how `resp` is read, and update(prior = ). Each against brms
# 2.23.0, run in dev/predread-me-compare.R and dev/predread-resp-brms-probe.R.
#
# Seen to fail on 0.25.0 (rellib-r16): the sd over draws of posterior_epred()
# on new data was 0.603 of brms's, the plug-in b0 + b x
# (dev/rel078-fix2-me-newdata.R); posterior_linpred() read a cs() model as
# the shared part alone; bayes_R2(resp = character(0)) stopped with "'data' must
# be of a vector type"; conditional_effects() of a multivariate model gave the
# first response; update(prior = ) dropped the first call's rows.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

rm_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

rd_data <- function() {
  set.seed(3)
  n <- 80
  g <- rep(1:8, each = 10)
  lat <- stats::rnorm(8)[g]
  d <- data.frame(x = stats::rnorm(n), g = factor(g), sg = 0.3, sx = 0.4)
  d$xg <- lat + stats::rnorm(8, 0, 0.3)[g]
  d$y <- 1 + 0.8 * d$x + rnorm(n, 0, 0.5)
  d$y2 <- -d$x + 0.5 * lat + rnorm(n, 0, 0.6)
  d$y3 <- 0.5 * d$x + rnorm(n)
  d
}

rd_cache <- new.env()
rd_fit <- function(key, ...) {
  if (is.null(rd_cache[[key]])) {
    rd_cache[[key]] <- rm_quiet(frm_sample(
      ..., data = rd_data(), family = gaussian(), chains = 1, iter = 300,
      seed = 3, refresh = 0))
  }
  rd_cache[[key]]
}

test_that("me() on new data draws a new latent value per draw", {
  ds <- rd_fit("me", bf(y ~ me(x, sx)))
  nd <- data.frame(x = c(-1, 0, 1), sx = c(0.2, 0.4, 0.8))
  m <- as.matrix(ds)
  b0 <- m[, "b_Intercept"]
  b1 <- m[, "bsp_mexsx"]
  set.seed(1)
  rp <- 40L
  e <- do.call(rbind, replicate(rp, posterior_epred(ds, newdata = nd),
                                simplify = FALSE))
  # the same draws, b0 + b (x + sdx z), with z fresh, brms's construction
  set.seed(2)
  want <- do.call(rbind, replicate(rp, {
    sapply(1:3, function(i) {
      b0 + b1 * (nd$x[i] + nd$sx[i] * stats::rnorm(length(b0)))
    })
  }, simplify = FALSE))
  n <- nrow(want)
  ratio <- apply(e, 2, stats::sd) / apply(want, 2, stats::sd)
  expect_lt(max(abs(ratio - 1)), 5 / sqrt(2 * n))
  zm <- (colMeans(e) - colMeans(want)) / (apply(want, 2, stats::sd) / sqrt(n))
  expect_lt(max(abs(zm)), 5)
  # it moves: the plug-in's spread is smaller on every row
  plug <- sapply(1:3, function(i) b0 + b1 * nd$x[i])
  expect_true(all(apply(e, 2, stats::sd) > 1.2 * apply(plug, 2, stats::sd)))
  # the absent case: on the fitted rows the draws' own latent values are read,
  # draw for draw
  expect_identical(posterior_epred(ds, draw_ids = 1:5),
                   posterior_epred(ds, draw_ids = 1:5))
  # posterior_predict() and posterior_linpred() read it too
  l <- posterior_linpred(ds, newdata = nd)
  expect_false(isTRUE(all.equal(l, posterior_linpred(ds, newdata = nd))))
  # a missing sdx is refused, as brms refuses it
  expect_error(posterior_epred(ds, newdata = nd["x"]), "measurement SD")
  expect_error(posterior_predict(ds, newdata = transform(nd, sx = 0)),
               "should be positive")
})

test_that("a new row with no x or sdx reads NA on draws, the others answer", {
  # brms reads NaN and warns; a refusal stopped the rows that have values
  ds <- rd_fit("me", bf(y ~ me(x, sx)))
  nd <- data.frame(x = c(-1, 0, 1), sx = c(0.2, NA, 0.8))
  for (f in list(posterior_epred, posterior_linpred)) {
    n <- 0L
    e <- withCallingHandlers(f(ds, newdata = nd), warning = function(w) {
      if (inherits(w, "frmtmb_me_missing")) n <<- n + 1L
      invokeRestart("muffleWarning")
    })
    # once in the call, whatever the number of draws
    expect_identical(n, 1L)
    expect_true(all(is.na(e[, 2L])))
    expect_false(anyNA(e[, -2L]))
  }
})

test_that("me() on new data keeps working under laplace = TRUE", {
  # the guard compares two fills under one stream, so a drawn latent value
  # does not read as a random effect; the fitted rows are refused as before
  ds <- rm_quiet(frm_sample(bf(y ~ me(x, sx)), data = rd_data(),
                            family = gaussian(), chains = 1, iter = 200,
                            seed = 3, refresh = 0, laplace = TRUE))
  nd <- data.frame(x = c(-1, 1), sx = 0.4)
  expect_identical(dim(posterior_epred(ds, newdata = nd))[2L], 2L)
  expect_error(posterior_epred(ds), "integrates them out")
})

test_that("me(gr = ) rows of one level share a draw", {
  ds <- rd_fit("gr", bf(y2 ~ me(xg, sg, gr = g)))
  nd <- data.frame(xg = c(0.5, 0.5, -1), sg = 0.3,
                   g = factor(c("1", "1", "2"), levels = as.character(1:8)))
  set.seed(1)
  e <- posterior_epred(ds, newdata = nd)
  expect_identical(e[, 1L], e[, 2L])
  expect_false(isTRUE(all.equal(e[, 1L], e[, 3L])))
  bad <- transform(nd, xg = c(0.5, 9, -1))
  expect_error(posterior_epred(ds, newdata = bad), "unique for each group")
})

test_that("posterior_linpred() of a cs() model has brms's threshold layers", {
  set.seed(4)
  d <- data.frame(x = stats::rnorm(150), z = stats::rnorm(150))
  d$o <- cut(d$x - 0.5 * d$z + stats::rlogis(150),
             c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
  ds <- rm_quiet(frm_sample(bf(o ~ x + cs(z)), data = d, family = sratio(),
                            chains = 1, iter = 300, seed = 3, refresh = 0))
  nd <- data.frame(x = c(-1, 1), z = c(0.5, -1))
  l <- posterior_linpred(ds, newdata = nd)
  # draws x rows x thresholds, every margin unnamed (brms 2.23.0,
  # dev/predread-cs-layout-brms.R)
  expect_identical(dim(l)[2:3], c(2L, 2L))
  expect_null(dimnames(l))
  m <- as.matrix(ds)
  want <- m[, "b_x"] * nd$x[1L] + m[, "bcs_z[1]"] * nd$z[1L]
  expect_equal(unname(l[, 1L, 1L]), unname(want))
  want2 <- m[, "b_x"] * nd$x[2L] + m[, "bcs_z[2]"] * nd$z[2L]
  expect_equal(unname(l[, 2L, 2L]), unname(want2))
  # dpar = "mu" is the shared part, and so is transform = TRUE's
  sh <- posterior_linpred(ds, newdata = nd, dpar = "mu")
  expect_identical(dim(sh), c(nrow(m), 2L))
  expect_equal(unname(sh[, 1L]), unname(m[, "b_x"] * nd$x[1L]))
  tr <- posterior_linpred(ds, newdata = nd, transform = TRUE)
  expect_identical(dim(tr), c(nrow(m), 2L))
  # the predictive draws still read the shared part beside .cs
  expect_identical(dim(posterior_predict(ds, newdata = nd)), c(nrow(m), 2L))
  ep <- posterior_epred(ds, newdata = nd)
  expect_identical(dim(ep)[3L], 3L)
  thr <- posterior_linpred(ds, newdata = nd, incl_thres = TRUE)
  tau <- m[, "b_Intercept[1]"]
  expect_equal(unname(thr[, 1L, 1L]), unname(tau - l[, 1L, 1L]))
})

test_that("an empty resp is every response in bayes_R2() and loo_R2()", {
  ds <- rd_fit("mv", mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE))
  expect_identical(rownames(bayes_R2(ds, resp = character(0))),
                   c("R2y", "R2y2"))
  expect_equal(bayes_R2(ds, resp = character(0)), bayes_R2(ds))
  skip_if_not_installed("loo")
  lr <- rm_quiet(loo_R2(ds, resp = character(0)))
  expect_identical(rownames(lr), c("R2y", "R2y2"))
  # a name given twice stays refused here; brms answers two equal rows
  expect_error(bayes_R2(ds, resp = c("y", "y")), "twice")
})

test_that("conditional_effects() gives one panel set per response", {
  ds <- rd_fit("mv", mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE))
  nm <- function(...) {
    names(rm_quiet(conditional_effects(ds, "x", ..., resolution = 3L)))
  }
  # as brms 2.23.0 names them
  expect_identical(nm(), c("y.y_x", "y2.y2_x"))
  expect_identical(nm(resp = character(0)), c("y.y_x", "y2.y2_x"))
  expect_identical(nm(resp = "y2"), "y2.y2_x")
  expect_identical(nm(resp = c("y2", "y")), c("y2.y2_x", "y.y_x"))
  expect_identical(nm(resp = c("y", "y2"), method = "posterior_predict"),
                   c("y.y_x", "y2.y2_x"))
  expect_error(nm(resp = c("y", "y")), "twice")
  expect_error(nm(resp = "nosuch"), "resp")
  ce <- rm_quiet(conditional_effects(ds, "x", resolution = 3L))
  expect_identical(attr(ce[["y2.y2_x"]], "response"), "y2")
  # the panel is its response's own
  ep <- posterior_epred(ds, newdata = data.frame(x = ce[["y2.y2_x"]]$x),
                        resp = "y2")
  expect_equal(ce[["y2.y2_x"]]$estimate__,
               unname(apply(ep, 2, stats::median)))
})

test_that("a univariate model does not read resp", {
  ds <- rd_fit("uni", bf(y ~ x))
  nd <- data.frame(x = c(-1, 1))
  want <- posterior_epred(ds, newdata = nd)
  wlin <- posterior_linpred(ds, newdata = nd)
  wfit <- fitted(ds, newdata = nd)
  wr2 <- bayes_R2(ds)
  # as brms 2.23.0, each run: an empty vector, a name given twice, an unknown
  # name, NA and a number all read as no resp
  for (r in list(character(0), c("y", "y"), "nosuch", NA_character_, 1L)) {
    expect_identical(posterior_epred(ds, newdata = nd, resp = r), want)
    expect_identical(posterior_linpred(ds, newdata = nd, resp = r), wlin)
    expect_identical(fitted(ds, newdata = nd, resp = r), wfit)
    expect_identical(dim(posterior_predict(ds, newdata = nd, resp = r)),
                     dim(want))
    expect_identical(dim(log_lik(ds, resp = r)), c(nrow(want), 80L))
    expect_identical(bayes_R2(ds, resp = r), wr2)
    expect_identical(names(rm_quiet(conditional_effects(
      ds, "x", resp = r, resolution = 3L))), "x")
  }
})

test_that("pp_mixture() reads one named response of a multivariate model", {
  d <- rd_data()
  ds <- rm_quiet(frm_sample(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE),
                            data = d, family = mixture(gaussian(), gaussian()),
                            chains = 1, iter = 200, seed = 3, refresh = 0,
                            init = 0))
  for (r in list(NULL, character(0), c("y", "y2"))) {
    expect_error(pp_mixture(ds, resp = r), "single variable name|twice")
  }
  expect_error(pp_mixture(ds, resp = c("y", "y")), "twice")
  expect_error(pp_mixture(ds, resp = "nosuch"), "resp")
  p1 <- pp_mixture(ds, resp = "y", draw_ids = 1:3)
  p2 <- pp_mixture(ds, resp = "y2", draw_ids = 1:3)
  expect_identical(dim(p1), c(nrow(d), 4L, 2L))
  expect_false(isTRUE(all.equal(p1, p2)))
})

test_that("update(prior = ) merges with the first call's prior, as brms", {
  ds <- rd_fit("pr", bf(y ~ x), prior = c(
    set_prior("normal(9, 0.01)", class = "b", coef = "x"),
    set_prior("normal(0, 10)", class = "Intercept")))
  up <- rm_quiet(update(ds, prior = set_prior("exponential(1)",
                                              class = "sigma")))
  ps <- as.data.frame(prior_summary(up))
  have <- function(cls, cf = "") {
    ps$prior[ps$class == cls & ps$coef == cf]
  }
  # the first call's rows are still in, and the new one is there
  expect_identical(have("b", "x"), "normal(9, 0.01)")
  expect_identical(have("Intercept"), "normal(0, 10)")
  expect_identical(have("sigma"), "exponential(1)")
  expect_lt(abs(mean(as.matrix(up)[, "b_x"]) - 9), 0.1)
  # a new row takes its slot
  up2 <- rm_quiet(update(ds, prior = set_prior("normal(-9, 0.01)", class = "b",
                                               coef = "x")))
  expect_lt(abs(mean(as.matrix(up2)[, "b_x"]) + 9), 0.1)
  ps2 <- as.data.frame(prior_summary(up2))
  expect_identical(ps2$prior[ps2$class == "Intercept" & ps2$coef == ""],
                   "normal(0, 10)")
  # a prior with no rows keeps them
  up3 <- rm_quiet(update(ds, prior = NULL, iter = 200))
  expect_lt(abs(mean(as.matrix(up3)[, "b_x"]) - 9), 0.1)
  # the absent case: another argument changes nothing about the prior
  up4 <- rm_quiet(update(ds, iter = 200))
  expect_lt(abs(mean(as.matrix(up4)[, "b_x"]) - 9), 0.1)
})

test_that("a multivariate model's cs() layers have brms's names", {
  set.seed(4)
  d <- data.frame(x = stats::rnorm(150), z = stats::rnorm(150))
  d$o <- cut(d$x - 0.5 * d$z + stats::rlogis(150),
             c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
  d$y <- d$x + stats::rnorm(150)
  ds <- rm_quiet(frm_sample(
    mvbf(bf(y ~ x), bf(o ~ x + cs(z), family = sratio()), rescor = FALSE),
    data = d, chains = 1, iter = 200, seed = 3, refresh = 0))
  nd <- data.frame(x = c(-1, 1), z = c(0.5, -1), y = 0)
  # brms 2.23.0 on the same call (dev/predread-sample-cells.R)
  l <- posterior_linpred(ds, newdata = nd)
  expect_identical(dimnames(l)[[3L]], c("y", "o1", "o2"))
  expect_identical(dim(posterior_linpred(ds, newdata = nd, resp = "o"))[2:3],
                   c(2L, 2L))
  fi <- fitted(ds, newdata = nd, scale = "linear")
  expect_identical(dimnames(fi)[[3L]], c("y", "eta1", "eta2"))
  # a univariate summary too
  du <- rm_quiet(frm_sample(bf(o ~ x + cs(z)), data = d, family = sratio(),
                            chains = 1, iter = 200, seed = 3, refresh = 0))
  expect_identical(dimnames(fitted(du, newdata = nd, scale = "linear"))[[3L]],
                   c("eta1", "eta2"))
})

test_that("responses of one call read one latent of a me() term per draw", {
  # a stated divergence: brms's own Stan program has one Xme_1, and its draws on
  # new data contradict it (dev/review-predread-mejoint.R)
  set.seed(7701)
  n <- 80
  lat <- stats::rnorm(n)
  d <- data.frame(xo = lat + stats::rnorm(n, 0, 0.4), sx = 0.4)
  d$y1 <- 1.5 * lat + stats::rnorm(n, 0, 0.4)
  d$y2 <- 1.2 * lat + stats::rnorm(n, 0, 0.4)
  ds <- rm_quiet(frm_sample(
    mvbf(bf(y1 ~ me(xo, sx)), bf(y2 ~ me(xo, sx)), rescor = FALSE), data = d,
    family = gaussian(), chains = 1, iter = 300, seed = 3, refresh = 0))
  nd <- data.frame(xo = c(0, 1), sx = c(1, 2))
  set.seed(1)
  e <- posterior_epred(ds, newdata = nd)
  expect_gt(stats::cor(e[, 2L, 1L], e[, 2L, 2L]), 0.97)
  pp <- posterior_predict(ds, newdata = nd)
  expect_gt(stats::cor(pp[, 2L, 1L], pp[, 2L, 2L]), 0.8)
  l <- posterior_linpred(ds, newdata = nd)
  expect_gt(stats::cor(l[, 2L, 1L], l[, 2L, 2L]), 0.97)
  # a call apart from it is another draw
  e2 <- posterior_epred(ds, newdata = nd, resp = "y2")
  expect_lt(abs(stats::cor(e[, 2L, 1L], e2[, 2L])), 0.3)
  # sigma and mu read one latent too
  ds2 <- rm_quiet(frm_sample(bf(y1 ~ me(xo, sx), sigma ~ me(xo, sx)), data = d,
                             family = gaussian(), chains = 1, iter = 300,
                             seed = 3, refresh = 0))
  sh <- frmtmb.sample:::draws_fit_at(ds2, 1L)
  set.seed(2)
  m1 <- frm_linpred(sh, newdata = nd, dpar = "mu")
  s1 <- frm_linpred(sh, newdata = nd, dpar = "sigma")
  ms <- as.matrix(ds2)[1L, ]
  z_m <- (m1 - ms[["b_Intercept"]]) / ms[["bsp_mexosx"]]
  z_s <- (s1 - ms[["b_sigma_Intercept"]]) / ms[["bsp_sigma_mexosx"]]
  expect_equal(unname(z_m), unname(z_s))
})

test_that("a display of me(gr = ) on draws holds the group at NA, and answers", {
  ds <- rd_fit("gr", bf(y2 ~ me(xg, sg, gr = g)))
  ce <- rm_quiet(conditional_effects(ds, "xg", resolution = 3L))[[1L]]
  expect_true(all(is.na(ce$g)))
  expect_true(all(is.finite(ce$estimate__)))
  ce2 <- rm_quiet(conditional_effects(ds, "xg", resolution = 3L,
                                      method = "posterior_predict"))[[1L]]
  expect_true(all(is.finite(ce2$estimate__)))
})

test_that("posterior_predict(newdata) answers under cov = TRUE", {
  set.seed(9)
  d <- data.frame(t = 1:60, x = stats::rnorm(60))
  d$y <- 1 + 0.5 * d$x + as.numeric(stats::arima.sim(list(ar = 0.5), 60))
  ds <- rm_quiet(frm_sample(bf(y ~ x + ar(time = t, cov = TRUE)), data = d,
                            family = gaussian(), chains = 1, iter = 300,
                            seed = 3, refresh = 0))
  nd <- data.frame(x = c(-1, 0, 1), t = 5:7)
  p <- posterior_predict(ds, newdata = nd)
  expect_identical(ncol(p), 3L)
  expect_gt(stats::cor(p[, 1L], p[, 2L]), 0.2)
})

test_that("panels on draws are in brms's order", {
  ds <- rd_fit("mv3", mvbf(bf(y ~ x), bf(y2 ~ x + y3), rescor = FALSE))
  ce <- rm_quiet(conditional_effects(ds, resolution = 3L))
  expect_named(ce, c("y.y_x", "y2.y2_x", "y.y_y3", "y2.y2_y3"))
})
