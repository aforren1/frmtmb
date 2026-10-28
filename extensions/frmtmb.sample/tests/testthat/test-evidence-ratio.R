# The Savage-Dickey evidence ratio and the posterior probability beside
# it, on the frmtmb_draws hypothesis() method.
#
# The identities are asserted against the construction rather than
# against remembered numbers: a directional ratio IS the posterior odds
# of the claim, and a point ratio IS brms's evidence_ratio(), the kernel
# density of the posterior draws at the point over the kernel density of
# the prior draws there. A remembered number would only say that this
# machine's chain has not moved. The agreement with brms's own numbers
# on one model is measured in dev/priordraws-evidence.R.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

er_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      set.seed(21)
      n <- 150
      dd <- data.frame(x = stats::rnorm(n), w = stats::rnorm(n),
                       g = factor(rep(1:10, 15)))
      dd$y <- stats::rnorm(n, 0.5 + 0.6 * dd$w +
                             stats::rnorm(10, 0, 0.5)[dd$g], 1)
      fit <- frm(bf(y ~ x + w + (1 | g)), family = gaussian(), data = dd)
      pr <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 1200, warmup = 400,
                   refresh = 0, seed = 7, sample_prior = "yes",
                   prior = set_prior("normal(0, 1)", class = "b"))))
      # the same posterior without prior draws
      np <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 1200, warmup = 400,
                   refresh = 0, seed = 7,
                   prior = set_prior("normal(0, 1)", class = "b"))))
      flat <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 600, refresh = 0, seed = 7,
                   sample_prior = "yes")))
      cache <<- list(fit = fit, pr = pr, np = np, flat = flat)
    }
    cache
  }
})

test_that("the frame carries brms's Evid.Ratio and Post.Prob", {
  cs <- er_case()
  h <- hypothesis(cs$pr, c("x = 0", "x > 0"))$hypothesis
  expect_true(all(c("Evid.Ratio", "Post.Prob") %in% names(h)))
  expect_true(all(is.finite(h$Evid.Ratio)))
  expect_equal(h$Post.Prob, h$Evid.Ratio / (1 + h$Evid.Ratio))
  # printing must survive the character columns: print() rounds the
  # numeric ones only, so a non-numeric one must not abort it
  expect_output(print(hypothesis(cs$pr, c("x = 0", "x > 0"))), "Evid.Ratio")
})

test_that("a directional ratio is the posterior odds of the claim", {
  cs <- er_case()
  ho <- hypothesis(cs$pr, c("x > 0", "x < 0"))
  h <- ho$hypothesis
  dr <- as.matrix(ho$samples)
  gt <- sum(dr[, 1] > 0)
  lt <- sum(dr[, 2] < 0)
  n <- nrow(dr)
  expect_equal(h$Evid.Ratio[1], gt / (n - gt))
  expect_equal(h$Evid.Ratio[2], lt / (n - lt))
  # a directional ratio needs no prior, so it is reported on a model
  # whose slopes are flat, and on one without prior draws
  hf <- hypothesis(cs$flat, "x > 0")$hypothesis
  expect_true(is.finite(hf$Evid.Ratio[1]))
  expect_true(is.finite(hypothesis(cs$np, "x > 0")$hypothesis$Evid.Ratio))
})

test_that("a point ratio is the posterior density over the prior draws'", {
  cs <- er_case()
  h <- hypothesis(cs$pr, "x = 0")
  post <- h$samples$H1
  pri <- h$prior_samples$H1
  # the prior draws of b_x are the class-wide normal(0, 1) column
  expect_true(all(is.finite(pri)))
  expect_identical(sort(pri), sort(prior_draws(cs$pr)$b))
  expect_equal(h$hypothesis$Evid.Ratio[1],
               frmtmb.sample:::er_kde_at(post) /
                 frmtmb.sample:::er_kde_at(pri))
  # the Monte Carlo error rides along as an attribute, so the frame
  # keeps the shape a ported brms script indexes: the spread of the
  # numerator across four blocks, over the one denominator
  expect_length(attr(h, "evid_ratio_mcse"), 1L)
  nd <- length(post)
  ix <- split(seq_len(nd), rep(seq_len(4L), each = nd %/% 4L,
                               length.out = nd))
  per <- vapply(ix, function(i) {
    frmtmb.sample:::er_kde_at(post[i])
  }, numeric(1)) / frmtmb.sample:::er_kde_at(pri)
  expect_equal(attr(h, "evid_ratio_mcse")[[1]], stats::sd(per) / sqrt(4))
})

test_that("a hypothesis of several parameters reads each one's draws", {
  cs <- er_case()
  # b_x and b_w both read the class-wide column, each in its own random
  # order, as brms's prior_draws() does, so the prior of x - w is the
  # difference of two independent normal(0, 1) draws
  h <- hypothesis(cs$pr, "x - w = 0", seed = 3)
  pri <- h$prior_samples$H1
  expect_true(all(is.finite(pri)))
  ks <- stats::ks.test(pri, function(q) stats::pnorm(q, 0, sqrt(2)))
  expect_gt(ks$p.value, 0.001)
  expect_equal(h$hypothesis$Evid.Ratio[1],
               frmtmb.sample:::er_kde_at(h$samples$H1) /
                 frmtmb.sample:::er_kde_at(pri))
})

test_that("brms gives a point ratio for sigma and sd, so does this", {
  cs <- er_case()
  # brms reads prior_sigma and prior_sd_g, by the longest name that
  # starts each variable. Through frmtmb.sample 0.23.0 these were
  # refused with a warning
  hs <- hypothesis(cs$pr, c("sigma = 0.9", "sd_g__Intercept = 0.5"),
                   class = NULL)
  expect_true(all(is.finite(hs$hypothesis$Evid.Ratio)))
  expect_equal(sort(hs$prior_samples$H1 + 0.9),
               sort(prior_draws(cs$pr)$sigma))
  # a nonlinear function of a coefficient has prior draws too
  hq <- hypothesis(cs$pr, "x^2 = 0")
  expect_true(all(is.finite(hq$prior_samples$H1)))
})

test_that("without prior draws a point ratio is NA, with no message", {
  cs <- er_case()
  # brms's evidence_ratio() without prior_samples. Through frmtmb.sample
  # 0.23.0 the denominator was read off the prior specification instead
  expect_no_warning(h <- hypothesis(cs$np, c("x = 0", "x > 0")))
  expect_true(is.na(h$hypothesis$Evid.Ratio[1]))
  expect_true(is.na(h$hypothesis$Post.Prob[1]))
  expect_true(all(is.na(h$prior_samples$H1)))
  # and a slope left flat has no prior draws even with "yes"
  expect_no_warning(hf <- hypothesis(cs$flat, "x = 0"))
  expect_true(is.na(hf$hypothesis$Evid.Ratio[1]))
})

test_that("an overall intercept has no prior draws, as in brms", {
  cs <- er_case()
  # brms's ignore_prior(): the centered intercept's prior is not the
  # prior of b_Intercept, and brms warns and leaves it out
  h <- allow_warnings(hypothesis(cs$pr, c("Intercept = 0", "x = 0")),
                      "prior of an overall intercept",
                      require = "prior of an overall intercept")
  expect_true(is.na(h$hypothesis$Evid.Ratio[1]))
  expect_true(is.finite(h$hypothesis$Evid.Ratio[2]))
})

test_that("the kernel density estimate reaches a point outside the draws", {
  # density() evaluates on [min, max] only; a hypothesis whose point is
  # past the drawn range would otherwise be read off the end of the grid
  set.seed(3)
  v <- stats::rnorm(4000, 5, 1)
  d0 <- frmtmb.sample:::er_kde_at(v, point = 0)
  expect_true(is.finite(d0))
  expect_gte(d0, 0)
  # and at the centre it estimates the density that made the draws.
  # A gaussian kernel of bandwidth h on normal(5, 1) draws has the
  # expectation dnorm(0, 0, sqrt(1 + h^2)) at 5 exactly, and the
  # standard error sqrt(f R(K) / (n h)), R(K) = 1 / (2 sqrt(pi)), to
  # first order; the estimate is within 4 of them
  h <- stats::bw.nrd0(v)
  ex <- stats::dnorm(0, 0, sqrt(1 + h^2))
  se <- sqrt(stats::dnorm(0) / (2 * sqrt(pi)) / (length(v) * h))
  expect_lt(abs(frmtmb.sample:::er_kde_at(v, point = 5) - ex) / se, 4)
})

test_that("the MCSE has at least four blocks whatever the chain count", {
  cs <- er_case()
  # two chains give sd() of two numbers over sqrt(2): one degree of
  # freedom, which understated a measured 0.46 gap as 0.18
  expect_identical(nchains(cs$pr), 2L)
  n <- nrow(cs$pr$draws)
  expect_length(frmtmb.sample:::er_blocks(cs$pr, n), 4L)
  # and a draws object with more chains keeps its own boundaries
  expect_length(frmtmb.sample:::er_blocks(
    structure(list(stanfit = NULL, draws = matrix(0, 800, 1),
                   fit = cs$fit), class = "frmtmb_draws"), 800L), 4L)
})

## ---- the numerator is brms code, asserted not just claimed ----------

test_that("er_kde_at() is brms density_ratio(), bit for bit", {
  # No Stan and no fit: this compares two smoothers on a fixed vector,
  # which is the whole content of the claim that the two packages
  # numbers are comparable at all. It needs brms installed and nothing
  # else.
  skip_unless_brms()
  dr <- utils::getFromNamespace("density_ratio", "brms")
  set.seed(404)
  v <- stats::rnorm(1600, 0.12, 0.09)
  expect_identical(frmtmb.sample:::er_kde_at(v), dr(v, point = 0))
  # away from the drawn range, where the widening branch fires
  w <- stats::rnorm(2000, 5, 1)
  expect_identical(frmtmb.sample:::er_kde_at(w), dr(w, point = 0))
  # and at a point inside it
  expect_identical(frmtmb.sample:::er_kde_at(w, point = 5),
                   dr(w, point = 5))
  # the ratio of two of them is brms's density_ratio(x, y)
  expect_identical(frmtmb.sample:::er_kde_at(v) /
                     frmtmb.sample:::er_kde_at(w),
                   dr(v, w, point = 0))
})

test_that("plot() of a hypothesis overlays the prior draws, as brms", {
  cs <- er_case()
  h <- hypothesis(cs$pr, "x = 0")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grDevices::dev.control("enable")
  expect_false(plot(h, plot = FALSE)[[1]]$ignore_prior)
  # the page with the prior's density has one more drawing operation
  # than the page without it
  print(plot(h, plot = FALSE)[[1]])
  with_prior <- length(grDevices::recordPlot()[[1]])
  print(plot(h, ignore_prior = TRUE, plot = FALSE)[[1]])
  without <- length(grDevices::recordPlot()[[1]])
  expect_identical(with_prior - without, 1L)
  # draws without prior draws have nothing to overlay
  print(plot(hypothesis(cs$np, "x = 0"), plot = FALSE)[[1]])
  expect_identical(length(grDevices::recordPlot()[[1]]), without)
})
