## Phase 0 of dev/extension-gaps-plan.md: the sample row.
##
## A 2000-row GLMM with two crossed factors, 4 chains x 2000 draws, then
## posterior_epred() and loo() on the draws.
##
## What the row decides: how slow the per-draw R loops are, and whether
## Phase 4's item 4.6 caching is needed before anything else in this
## package. So the post-fit calls are timed one at a time, and
## posterior_epred() is timed at two draw counts, because a per-draw
## loop and a cached design differ in their SLOPE in the number of draws
## and not only in their level. A loop that is linear in draws with a
## large intercept is a different item from one that is linear with a
## large slope, and only the second is what 4.6 removes.
##
## "A 2000-row GLMM" does not name a family; binomial with a logit link
## is chosen here, because that is the model a sampler is reached for
## when the Laplace approximation is the thing in doubt, and because a
## gaussian response would make the posterior nearly normal and the
## sampler's cost unrepresentative.
##
## See dev/scale-findings.md for the numbers this produced.

sample_truth <- list(n = 2000L, ng1 = 40L, ng2 = 25L,
                     b0 = -0.3, b_x = 0.8, sd_g1 = 0.7, sd_g2 = 0.4,
                     chains = 4L, iter = 2000L)

sample_scale_data <- function(seed = 20260908L) {
  tr <- sample_truth
  set.seed(seed)
  n <- if (scale_small()) 200L else tr$n
  d <- data.frame(x = stats::rnorm(n),
                  g1 = factor(rep_len(seq_len(tr$ng1), n)),
                  g2 = factor(sample.int(tr$ng2, n, replace = TRUE)))
  u1 <- stats::rnorm(tr$ng1, 0, tr$sd_g1)
  u2 <- stats::rnorm(tr$ng2, 0, tr$sd_g2)
  eta <- tr$b0 + tr$b_x * d$x + u1[as.integer(d$g1)] +
    u2[as.integer(d$g2)]
  d$y <- stats::rbinom(n, 1L, stats::plogis(eta))
  d
}

test_that("the sample scale row samples and reports its cost", {
  skip_unless_scale()
  scale_row_on("sample")
  skip_sampler()
  tr <- sample_truth
  d <- sample_scale_data()
  chains <- if (scale_small()) 2L else tr$chains
  iter <- if (scale_small()) 300L else tr$iter
  form <- bf(y ~ x + (1 | g1) + (1 | g2))
  scale_mem_reset()

  fit <- NULL
  t_ml <- scale_elapsed(fit <- frm(form, family = stats::binomial(),
                                   data = d, se = TRUE))
  # The four chains run SERIALLY. frm_sample() passes tmbstan's `cores`
  # through and rstan defaults it to getOption("mc.cores", 1), which is
  # 1 in the process this tier runs in, so `sample_s` is the cost of
  # four chains one after another and not of four at once.
  ds <- NULL
  t_sample <- scale_elapsed(
    ds <- suppressWarnings(suppressMessages(
      frm_sample(fit, chains = chains, iter = iter, refresh = 0,
                 seed = 1))))
  nd <- ndraws(ds)

  # INTERLEAVED and replicated, because the number this row exists to
  # produce is a RATIO between the first two arms, and a ratio of two
  # single timings is not a measurement.
  tenth <- max(1L, nd %/% 10L)
  pf <- scale_interleave(list(
    epred = function() posterior_epred(ds),
    epred_tenth = function() posterior_epred(ds, ndraws = tenth),
    predict = function() posterior_predict(ds),
    ranef = function() ranef(ds),
    log_lik = function() log_lik(ds),
    loo = function() suppressWarnings(loo(ds))))
  ep <- posterior_epred(ds)
  ll <- log_lik(ds)
  lo <- suppressWarnings(loo(ds))
  mem <- scale_mem_peak_mb()

  ps <- posterior_summary(ds)
  # the slope's label is whatever this package spells it, found rather
  # than assumed, and recorded so the row says which column it read
  b_x <- rownames(ps)[grepl("(^|_)x$", rownames(ps))][1L]
  scale_record(
    "sample", rows = nrow(d), draws = nd,
    chains = chains, iter = iter,
    ml_fit_s = t_ml, sample_s = t_sample,
    rounds = pf$rounds,
    epred_first_s = pf$first[["epred"]],
    epred_s = pf$seconds[["epred"]],
    epred_tenth_s = pf$seconds[["epred_tenth"]],
    epred_draw_ratio = nd / tenth,
    epred_time_ratio = pf$seconds[["epred"]] /
      pf$seconds[["epred_tenth"]],
    epred_spread = pf$spread[["epred"]],
    epred_tenth_spread = pf$spread[["epred_tenth"]],
    epred_per_draw_ms = 1000 * pf$seconds[["epred"]] / nd,
    predict_s = pf$seconds[["predict"]],
    ranef_s = pf$seconds[["ranef"]],
    log_lik_s = pf$seconds[["log_lik"]],
    loo_s = pf$seconds[["loo"]], mem_mb = mem,
    epred_dim = paste(dim(ep), collapse = "x"),
    log_lik_dim = paste(dim(ll), collapse = "x"),
    b_x_name = b_x, b_x_mean = ps[b_x, "Estimate"],
    b_x_lo = ps[b_x, 3L], b_x_hi = ps[b_x, 4L], b_x_true = tr$b_x,
    # every summarized parameter by name, so that a reader can see
    # which column the slope was read from and what the variance
    # components came out at
    par_summary = paste(rownames(ps),
                        formatC(ps[, "Estimate"], digits = 4,
                                format = "g"),
                        sep = "=", collapse = ";"),
    sd_g1_true = tr$sd_g1, sd_g2_true = tr$sd_g2,
    looic = lo$estimates["looic", "Estimate"],
    n_pareto_bad = sum(loo::pareto_k_values(lo) > 0.7),
    diag = scale_diag(fit))

  expect_true(all(is.finite(ep[1L, ])))
  expect_identical(dim(ll), c(nd, nrow(d)))
  # the slope is what the design is powered for; the assertion is the
  # posterior's own 95 percent interval covering the simulator's truth
  expect_true(ps[b_x, 3L] <= tr$b_x && ps[b_x, 4L] >= tr$b_x)
})
