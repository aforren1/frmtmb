# A partial re_formula on DRAWS keeps the terms it names, as brms's
# does. The draws methods evaluate each posterior draw through
# frmtmb::frm_linpred(), which read any formula as "keep every term"
# before frmtmb's development version; the fix is in frmtmb, and this
# file pins that it reaches every method on draws. Each draw carries
# its own sampled b, so the group-effect uncertainty at a known level is
# in these draws already and nothing is added here.

reunc_draws <- local({
  cache <- new.env(parent = emptyenv())
  function() {
    if (!is.null(cache$ds)) return(cache$ds)
    set.seed(3)
    G <- 10
    m <- 6
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    h = factor(rep(1:4, length.out = G * m)),
                    x = stats::rnorm(G * m))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(G, 0, 0.8)[d$g] +
      stats::rnorm(G, 0, 0.3)[d$g] * d$x + stats::rnorm(4, 0, 0.6)[d$h] +
      stats::rnorm(G * m)
    fit <- frm(bf(y ~ x + (1 + x | g) + (1 | h)), data = d)
    cache$d <- d
    cache$ds <- suppressWarnings(frm_sample(fit, chains = 1, iter = 160,
                                            warmup = 100, seed = 3,
                                            refresh = 0))
    cache$ds
  }
})

# one draw's expected response, built by hand from its own parameters
reunc_hand <- function(ds, k, ig, sl, ih) {
  sh <- frmtmb.sample:::draws_fit_at(ds, k)
  d <- ds$fit$frame[["data_frame"]]
  fe <- sh$estimates[["beta"]]
  B <- sh$estimates[["b"]]
  bk <- sh$frame[["re_blocks"]]
  bg <- t(matrix(B[bk[[1]][["b_idx"]]], 2))
  bh <- B[bk[[2]][["b_idx"]]]
  fe[1] + fe[2] * d$x + ig * bg[as.integer(d$g), 1] +
    sl * bg[as.integer(d$g), 2] * d$x + ih * bh[as.integer(d$h)]
}

test_that("posterior_epred() on draws keeps the named terms", {
  skip_on_cran()
  skip_sampler()
  ds <- reunc_draws()
  cases <- list(list(f = ~(1 | g), k = c(1, 0, 0)),
                list(f = ~(0 + x | g), k = c(0, 1, 0)),
                list(f = ~(1 | h), k = c(0, 0, 1)),
                list(f = ~(1 | g) + (1 | h), k = c(1, 0, 1)))
  for (cs in cases) {
    ep <- posterior_epred(ds, re_formula = cs$f)
    for (k in c(1L, nrow(ep))) {
      expect_equal(unname(ep[k, ]),
                   unname(reunc_hand(ds, k, cs$k[1], cs$k[2], cs$k[3])),
                   info = deparse1(cs$f))
    }
    lin <- posterior_linpred(ds, re_formula = cs$f)
    expect_equal(lin, ep, info = deparse1(cs$f))
    fi <- fitted(ds, re_formula = cs$f)
    expect_equal(unname(fi[, "Estimate"]), unname(colMeans(ep)))
  }
})

test_that("posterior_predict() and predict() on draws keep them too", {
  skip_on_cran()
  skip_sampler()
  ds <- reunc_draws()
  set.seed(5)
  pp <- posterior_predict(ds, re_formula = ~(1 | h))
  set.seed(5)
  pr <- predict(ds, re_formula = ~(1 | h))
  expect_equal(unname(pr[, "Estimate"]), unname(colMeans(pp)))
  # a gaussian draw is centred on its draw's expected response, so the
  # residual of the draws against the kept-terms epred has mean zero
  ep <- posterior_epred(ds, re_formula = ~(1 | h))
  res <- pp - ep
  z <- colMeans(res) / (apply(res, 2, stats::sd) / sqrt(nrow(res)))
  expect_lt(mean(abs(z) > 4), 0.05)
  # against the full prediction the same residual is off by the dropped
  # g terms, which is the defect this file pins
  full <- posterior_epred(ds)
  z_full <- colMeans(pp - full) /
    (apply(pp - full, 2, stats::sd) / sqrt(nrow(pp)))
  expect_gt(mean(abs(z_full) > 4), 0.05)
})

test_that("a term the fit does not have is refused on draws", {
  skip_on_cran()
  skip_sampler()
  ds <- reunc_draws()
  expect_error(posterior_epred(ds, re_formula = ~(1 | nosuch)),
               class = "frmtmb_error", regexp = "nosuch")
  expect_error(posterior_predict(ds, re_formula = ~(1 | nosuch)),
               class = "frmtmb_error", regexp = "nosuch")
  expect_error(fitted(ds, re_formula = ~(1 | nosuch)),
               class = "frmtmb_error", regexp = "nosuch")
})
