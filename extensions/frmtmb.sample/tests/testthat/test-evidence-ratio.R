# The Savage-Dickey evidence ratio and the posterior probability beside
# it, on the frmtmb_draws hypothesis() method.
#
# The identities are asserted against the construction rather than
# against remembered numbers: a directional ratio IS the posterior odds
# of the claim, and a point ratio IS the kernel density at the point
# over the prior density there. A remembered number would only say that
# this machine's chain has not moved.

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
                   refresh = 0, seed = 7,
                   prior = set_prior("normal(0, 1)", class = "b"))))
      flat <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 600, refresh = 0, seed = 7)))
      cache <<- list(fit = fit, pr = pr, flat = flat)
    }
    cache
  }
})

test_that("the frame carries evid_ratio and post_prob", {
  cs <- er_case()
  h <- hypothesis(cs$pr, c("x = 0", "x > 0"))
  expect_true(all(c("evid_ratio", "post_prob") %in% names(h)))
  expect_true(all(is.finite(h$evid_ratio)))
  expect_equal(h$post_prob, h$evid_ratio / (1 + h$evid_ratio))
  # printing must survive the new columns: print() rounds every column
  # after the first, so a non-numeric one would abort there
  expect_output(print(h), "evid_ratio")
})

test_that("a directional ratio is the posterior odds of the claim", {
  cs <- er_case()
  h <- hypothesis(cs$pr, c("x > 0", "x < 0"))
  dr <- attr(h, "draws")
  gt <- sum(dr[, 1] > 0)
  lt <- sum(dr[, 2] < 0)
  n <- nrow(dr)
  expect_equal(h$evid_ratio[1], gt / (n - gt))
  expect_equal(h$evid_ratio[2], lt / (n - lt))
  # a directional ratio needs no prior, so it is reported on a model
  # whose slopes are flat
  hf <- hypothesis(cs$flat, "x > 0")
  expect_true(is.finite(hf$evid_ratio[1]))
})

test_that("a point ratio is the density at zero over the prior there", {
  cs <- er_case()
  h <- hypothesis(cs$pr, "x = 0")
  dr <- attr(h, "draws")[, 1]
  expect_equal(h$evid_ratio[1],
               frmtmb.sample:::er_kde_at(dr) / stats::dnorm(0, 0, 1))
  # the Monte Carlo error rides along as an attribute, so the frame
  # keeps the shape a ported brms script indexes. Pinned to its
  # ESTIMATOR rather than to "smaller than the value", which any MCSE
  # under 400% of the ratio satisfies.
  expect_length(attr(h, "evid_ratio_mcse"), 1L)
  nd <- length(dr)
  ix <- split(seq_len(nd), rep(seq_len(4L), each = nd %/% 4L,
                               length.out = nd))
  per <- vapply(ix, function(i) {
    frmtmb.sample:::er_kde_at(dr[i]) / stats::dnorm(0, 0, 1)
  }, numeric(1))
  expect_equal(attr(h, "evid_ratio_mcse")[[1]],
               stats::sd(per) / sqrt(4))
})

test_that("an affine contrast convolves its normal priors", {
  cs <- er_case()
  h <- hypothesis(cs$pr, "x - w = 0")
  # x - w under two independent normal(0, 1) priors is normal(0, sqrt2)
  expect_equal(h$evid_ratio[1],
               frmtmb.sample:::er_kde_at(attr(h, "draws")[, 1]) /
                 stats::dnorm(0, 0, sqrt(2)))
  # and a scaled single coefficient transforms by the slope
  h2 <- hypothesis(cs$pr, "2 * x = 0")
  expect_equal(h2$evid_ratio[1],
               frmtmb.sample:::er_kde_at(attr(h2, "draws")[, 1]) /
                 (stats::dnorm(0, 0, 1) / 2))
})

test_that("a point hypothesis with no proper prior refuses by name", {
  cs <- er_case()
  # frm_sample() leaves class "b" flat by default, exactly as brms does
  expect_warning(h <- hypothesis(cs$flat, "x = 0"),
                 "x has no proper prior")
  expect_true(is.na(h$evid_ratio[1]))
  expect_true(is.na(h$post_prob[1]))
})

test_that("what a Savage-Dickey denominator cannot be is named", {
  cs <- er_case()
  # a variance component: its prior sits on the log standard deviation,
  # and a point null at zero is on the boundary anyway
  expect_warning(hypothesis(cs$pr, "sd_g__Intercept = 0"),
                 "natural-scale summary")
  # sigma IS a coefficient in this model, reported on another scale;
  # telling the user it is "not a population-level coefficient" was
  # true of the name and false of the model
  expect_warning(hypothesis(cs$pr, "sigma = 0.9"),
                 "natural-scale summary")
  # a nonlinear function of the coefficients has no closed-form prior
  expect_warning(hypothesis(cs$pr, "x^2 = 0"), "not an affine function")
  # the Intercept prior is about the intercept at the predictor means,
  # which is not the coefficient the hypothesis names
  expect_warning(hypothesis(cs$pr, "Intercept = 0"),
                 "not written about that coefficient itself")
  # every refusal still returns the rest of the row
  h <- suppressWarnings(hypothesis(cs$pr, c("x^2 = 0", "x > 0")))
  expect_true(is.na(h$evid_ratio[1]))
  expect_true(is.finite(h$evid_ratio[2]))
})

test_that("the kernel density estimate reaches a point outside the draws", {
  # density() evaluates on [min, max] only; a hypothesis whose point is
  # past the drawn range would otherwise be read off the end of the grid
  set.seed(3)
  v <- stats::rnorm(4000, 5, 1)
  d0 <- frmtmb.sample:::er_kde_at(v, point = 0)
  expect_true(is.finite(d0))
  expect_gte(d0, 0)
  # and at the centre it is close to the density that made the draws
  expect_lt(abs(frmtmb.sample:::er_kde_at(v, point = 5) -
                  stats::dnorm(0)), 0.02)
})


## ---- what the affine probe has to refuse ----------------------------

test_that("a hypothesis that is not affine in BOTH directions refuses", {
  cs <- er_case()
  # abs() is the identity on every non-negative point, and the fitting
  # points (0 and the unit vectors) are all non-negative, so a
  # one-sided check passed it. The prior of |X| for X ~ N(0, 1) is a
  # half-normal with twice the density at zero, so the ratio came back
  # exactly 2x too large with no warning.
  expect_warning(h <- hypothesis(cs$pr, "abs(x) = 0"),
                 "not an affine function")
  expect_true(is.na(h$evid_ratio[1]))
  # the negative probe must not cost a legitimate affine hypothesis
  expect_no_warning(hypothesis(cs$pr, c("x = 0", "2 * x = 0",
                                        "x - w = 0", "-x = 0")))
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
})
