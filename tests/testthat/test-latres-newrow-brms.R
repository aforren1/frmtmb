# A seen level's Est.Error on a new row beside a latent residual, against brms 2.23.0's draws.
# It read 4.21 and 21.35 on the first build of 0.79.0, where brms reads 1.31 and 5.85
# (dev/reviews/2026-11-04-release.md, finding 1).

skip_unless_brms()
skip_on_cran()
if (!identical(Sys.getenv("FRMTMB_BRMS_FIT_TESTS"), "true")) {
  skip("set FRMTMB_BRMS_FIT_TESTS=true to run the fits against brms")
}

test_that("a seen level's Est.Error is brms's", {
  set.seed(401)
  ng <- 8
  nt <- 25
  d <- expand.grid(t = 1:nt, s = factor(1:ng))
  d$x <- stats::rnorm(nrow(d))
  re <- stats::rnorm(ng, 0, 0.8)
  e <- stats::ave(stats::rnorm(nrow(d)), d$s, FUN = function(z) {
    stats::filter(z, 0.5, "recursive")
  })
  d$y <- stats::rpois(nrow(d), exp(3 + 0.3 * d$x + re[d$s] + 0.2 * e))
  d$y2 <- stats::rpois(nrow(d), exp(3 + 0.3 * d$x + re[d$s] + 0.06 * e))
  nd <- data.frame(t = 26, s = factor(c(1, 2)), x = c(0, 1))
  form <- y2 ~ x + (1 | s) + ar(time = t, gr = s, cov = TRUE)
  fit <- suppressWarnings(frm(form, data = d, family = poisson()))
  ours <- fitted(fit, newdata = nd)
  b <- suppressWarnings(suppressMessages(brms::brm(
    form, data = d, family = poisson(), chains = 2, iter = 1500, warmup = 750, seed = 8,
    refresh = 0, silent = 2)))
  ep <- brms::posterior_epred(b, newdata = nd)
  ref <- apply(ep, 2L, stats::sd)
  expect_lt(max(abs(ours[, "Est.Error"] / ref - 1)), 0.15)
  # the interval, on the row where it is widest
  q <- apply(ep, 2L, stats::quantile, c(0.025, 0.975))
  expect_lt(max(abs(ours[, "Q97.5"] / q[2L, ] - 1)), 0.15)
})
