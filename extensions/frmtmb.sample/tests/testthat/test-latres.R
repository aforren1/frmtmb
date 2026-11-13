# frm_sample() on brms's latent residuals: a family with no residual y - mu.
# The ML route's identity with brms is core's test-latres-brms.R; here, that the
# sampling route samples that model, under brms's defaults and names.

lrs_data <- function() {
  set.seed(167)
  d <- expand.grid(t = 1:5, g = factor(1:20))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  e <- as.vector(replicate(20, stats::arima.sim(list(ar = 0.5), 5, sd = 0.7)))
  d$cnt <- stats::rpois(n, exp(0.8 + 0.4 * d$x + e))
  d
}

test_that("the default priors are brms's def_scale_prior on sderr", {
  d <- lrs_data()
  p <- default_prior(bf(cnt ~ x + ar(t, g, cov = TRUE)), family = poisson(),
                     data = d, route = "sample")
  expect_match(p$prior[p$class == "sderr"], "student_t(3, 0, 2.5)",
               fixed = TRUE)
  expect_equal(p$lb[p$class == "sderr"], 0)
  expect_equal(p$lb[p$class == "ar"], -1)
  expect_equal(p$ub[p$class == "ar"], 1)
  pc <- default_prior(bf(cnt ~ x + cosy(t, g)), family = poisson(),
                      data = d, route = "sample")
  expect_equal(pc$lb[pc$class == "cosy"], 0)
  expect_equal(pc$ub[pc$class == "cosy"], 1)
  pu <- default_prior(bf(cnt ~ x + unstr(t, g)), family = poisson(),
                      data = d, route = "sample")
  expect_match(pu$prior[pu$class == "cortime"], "lkj(1)", fixed = TRUE)
})

test_that("the draws carry brms's variables and read the latents", {
  skip_on_cran()
  d <- lrs_data()
  dr <- suppressWarnings(suppressMessages(
    frm_sample(bf(cnt ~ x + ar(t, g, cov = TRUE)), family = poisson(), data = d,
               chains = 1, iter = 300, warmup = 150, seed = 11, refresh = 0)))
  v <- variables(dr)
  expect_true(all(c("b_Intercept", "b_x", "ar[1]", "sderr") %in% v))
  dd <- posterior::as_draws_df(dr)
  expect_true(all(abs(dd[["ar[1]"]]) < 1))
  expect_true(all(dd[["sderr"]] > 0))
  # the expected counts read each draw's latents
  ep <- posterior_epred(dr)
  expect_equal(dim(ep), c(150L, nrow(d)))
  fit <- dr$fit
  bk <- frame_ac_terms(fit$frame)[[1L]]
  expect_true(isTRUE(bk$latent))
  ll <- log_lik(dr)
  expect_equal(dim(ll), c(150L, nrow(d)))
  # each row's value is the family's density at that draw's mean
  expect_lt(max(abs(ll[1L, ] - stats::dpois(d$cnt, ep[1L, ], log = TRUE))),
            1e-8 * max(1, max(abs(ll[1L, ]))))
  pp <- posterior_predict(dr)
  expect_equal(dim(pp), c(150L, nrow(d)))
  # new rows draw new latents, so the draws of a row's mean spread by more
  # than its fixed part does
  nd <- d[d$g == "1", ]
  nd$g <- factor("n")
  pn <- posterior_predict(dr, newdata = nd)
  expect_equal(dim(pn), c(150L, nrow(nd)))
  en <- posterior_epred(dr, newdata = nd)
  dd2 <- posterior::as_draws_df(dr)
  fixed <- exp(outer(dd2[["b_Intercept"]], rep(1, nrow(nd))) +
                 outer(dd2[["b_x"]], nd$x))
  expect_gt(stats::sd(log(en[, 1L]) - log(fixed[, 1L])), 0.1)
  # a grid drops the term, as brms's incl_autocor = FALSE does
  ce <- conditional_effects(dr, "x")
  expect_true(is.data.frame(ce[[1L]]))
})

test_that("frm_sample() fits the forms frm() refuses", {
  skip_on_cran()
  d <- lrs_data()
  d$pos <- stats::rgamma(nrow(d), 4, 4 / exp(0.5 + 0.2 * d$x))
  expect_error(frm(bf(pos ~ x + ar(t, g)), family = Gamma(link = "log"),
                   data = d), "frm_sample()", fixed = TRUE)
  dr <- suppressWarnings(suppressMessages(
    frm_sample(bf(pos ~ x + ar(t, g)), family = Gamma(link = "log"), data = d,
               chains = 1, iter = 200, warmup = 100, seed = 11, refresh = 0)))
  vars <- variables(dr)
  expect_true(all(c("b_Intercept", "b_x", "ar[1]", "sderr") %in% vars))
  # cov = FALSE's states are brms's zerr, and b[k] is left to brms's slopes
  expect_identical(sum(grepl("^zerr\\[", vars)), nrow(d))
  expect_false(any(grepl("^b\\[", vars)))
})

test_that("cov = FALSE refuses new rows on draws, in brms's words", {
  skip_on_cran()
  d <- lrs_data()
  dr <- suppressWarnings(suppressMessages(
    frm_sample(bf(cnt ~ x + ar(t, g)), family = poisson(), data = d,
               chains = 1, iter = 200, warmup = 100, seed = 11, refresh = 0)))
  nd <- d[d$g == "1", ]
  w <- "Cannot predict new latent residuals when using cov = FALSE"
  expect_error(posterior_epred(dr, newdata = nd), w)
  expect_error(posterior_predict(dr, newdata = nd), w)
})

test_that("mi() on a latent response beside a latent residual samples", {
  skip_on_cran()
  # frm() refuses the pair, and the sampler holds every value
  set.seed(12301)
  d <- expand.grid(t = 1:8, g = factor(1:10))
  d <- d[order(d$g, d$t), ]
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rgamma(nrow(d), 5, 5 / exp(0.3 * d$x))
  d$y[c(3L, 20L, 41L)] <- NA
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y | mi() ~ x + ar(t, g, cov = TRUE)),
               family = Gamma(link = "log"), data = d, chains = 1, iter = 300,
               refresh = 0, seed = 19)))
  m <- as.matrix(ds$draws)
  expect_true(all(paste0("Ymi[", 1:3, "]") %in% colnames(m)))
  expect_true(all(c("ar[1]", "sderr") %in% colnames(m)))
  expect_false(anyNA(posterior_epred(ds)))
})
