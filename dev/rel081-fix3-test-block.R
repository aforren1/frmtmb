
# A GEV row outside a draw's support: log_lik() of draws reads the exact
# density through row_lpdf(exact = TRUE), -Inf, where the fitted
# objective keeps a finite barrier (-2.1e82 and below at these cells,
# dev/rel081-log/loo-exact-r20.txt). The cells come from draws moved by
# hand: sampled draws put no row there (0 of 44,000 in 22 fits that read
# the barrier, dev/reviews/2026-11-10-release.md, finding 2). The fitted
# rows as newdata reach the same reader. loo(), waic() and psis() stop
# with one classed error; through 0.80.0 loo() and psis() stopped in the
# loo package and waic() gave an elpd of NaN.
gev_off_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
      set.seed(349)
      d <- data.frame(x = stats::rnorm(40))
      d$y <- rg(40, 1 + 0.3 * d$x, 1, 0.3)
      fit <- suppressWarnings(frm(y ~ x, family = gen_extreme_value(),
                                  data = d))
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 1, iter = 400, seed = 349, refresh = 0)))
      # three draws with xi above 0 get the mean that puts the smallest y
      # below their lower end: t = 1 + xi (y - mu) / sigma = -xi / sigma
      i <- which.min(d$y)
      dr <- ds$draws
      k <- which(dr[, "xi"] > 0.05)[1:3]
      for (j in k) {
        mu <- d$y[i] + dr[j, "sigma"] / dr[j, "xi"] + 1
        dr[j, "b_Intercept"] <- mu - dr[j, "b_x"] * d$x[i]
      }
      moved <- ds
      moved$draws <- dr
      cache <<- list(d = d, ds = ds, moved = moved, i = i, k = k)
    }
    cache
  }
})

test_that("log_lik() of draws is -Inf off the GEV support, with newdata too", {
  cs <- gev_off_case()
  ll <- log_lik(cs$moved)
  expect_false(anyNA(cs$k))
  expect_true(all(ll[cs$k, cs$i] == -Inf))
  # the other draws are the sampled ones, and finite everywhere
  expect_true(all(is.finite(ll[-cs$k, ])))
  expect_identical(ll[-cs$k, ], log_lik(cs$ds)[-cs$k, ])
  # brms's newdata = fit$data and the data fitted reach the same reader
  expect_identical(log_lik(cs$moved, newdata = cs$moved$data), ll)
  expect_identical(log_lik(cs$moved, newdata = cs$d), ll)
})

test_that("loo(), waic() and psis() stop on a log_lik() cell that is -Inf", {
  cs <- gev_off_case()
  skip_if_not_installed("loo")
  ll <- log_lik(cs$moved)
  ncol_bad <- sum(colSums(!is.finite(ll)) > 0L)
  counts <- sprintf("in %d cells, at %d of %d observations and 3 of %d draws",
                    sum(!is.finite(ll)), ncol_bad, ncol(ll), nrow(ll))
  calls <- list(
    loo = function(nd) loo(cs$moved, newdata = nd),
    waic = function(nd) waic(cs$moved, newdata = nd),
    psis = function(nd) psis(cs$moved, nd))
  for (nm in names(calls)) {
    for (nd in list(NULL, cs$moved$data)) {
      e <- expect_error(suppressWarnings(calls[[nm]](nd)),
                        class = "frmtmb_loglik_not_finite", info = nm)
      expect_match(conditionMessage(e), counts, fixed = TRUE, info = nm)
      expect_match(conditionMessage(e), "outside the support",
                   fixed = TRUE, info = nm)
    }
  }
  # the same draws without the moved ones: the guard's condition absent
  keep <- setdiff(seq_len(nrow(ll)), cs$k)
  expect_s3_class(suppressWarnings(waic(cs$moved, draw_ids = keep)), "waic")
  expect_s3_class(suppressWarnings(loo(cs$moved, draw_ids = keep)), "loo")
  expect_s3_class(suppressWarnings(loo(cs$ds)), "loo")
})
