# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 6: latent
# mixtures.
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/ParameterEstimation/LatentMixtures.
# Sampling statements become `target += ..._lpdf(...)`; the vague
# regularizers are dropped (see test-bcm-data-analysis.R); and the
# ports' hand-written `log_sum_exp(lp_parts)` becomes `log_mix`, which
# is the same arithmetic with the mixing weight in front of it.
#
# The frmtmb spelling of a per-row latent class is mixture(). Two
# things the book's models need and a plain mixture() call does not
# offer are done with an OFFSET-ONLY sub-formula:
#
#   mu2 ~ 0 + offset(zero)     pins a component's rate at the link's
#                              origin, which on a logit is one half
#   theta1 ~ 0 + offset(zero)  pins the mixing weight at one half,
#                              because theta is a multinomial logit
#                              with the last component as reference
#
# A formula with no columns and an offset is a linear predictor that is
# entirely known, which is exactly what "this parameter is fixed at a
# value the book chose" means.

# ---------------------------------------------------------------------
# Exams_1 and Exams_2: 15 people, 40 questions. Some were guessing (a
# rate of one half) and the rest knew something.
# ---------------------------------------------------------------------

bcm_exams_data <- function() {
  data.frame(k = c(21, 17, 21, 18, 22, 31, 31, 34, 34, 35, 35, 36, 39,
                   36, 35),
             n = 40L, zero = 0,
             id = factor(1:15))
}

bcm_exams1_code <- function() {
  paste(
    "data {",
    "  int<lower=1> p; array[p] int<lower=0> k; int<lower=1> n;",
    "}",
    "parameters {",
    "  real<lower=0.5, upper=1> phi;",
    "}",
    "model {",
    "  for (i in 1 : p)",
    "    target += log_mix(0.5, binomial_lpmf(k[i] | n, phi),",
    "                           binomial_lpmf(k[i] | n, 0.5));",
    "}",
    sep = "\n")
}

# Exams_2 puts a rate on each person in the second group. The original
# draws it from a normal truncated to (0, 1); frmtmb's random effects
# are normal on the LINEAR PREDICTOR, which here is the logit, so the
# group distribution is logit-normal. That is a different prior for the
# same model and the adaptation is stated rather than hidden; a
# truncated-normal random effect is not in this grammar.
bcm_exams2_code <- function() {
  paste(
    "data {",
    "  int<lower=1> p; array[p] int<lower=0> k; int<lower=1> n;",
    "}",
    "parameters {",
    "  real mu; real<lower=0> sigma; vector[p] philogit;",
    "}",
    "model {",
    "  target += normal_lpdf(philogit | mu, sigma);",
    "  for (i in 1 : p)",
    "    target += log_mix(0.5,",
    "                      binomial_lpmf(k[i] | n, inv_logit(philogit[i])),",
    "                      binomial_lpmf(k[i] | n, 0.5));",
    "}",
    sep = "\n")
}

bcm_exams1_formula <- function() {
  bf(k | trials(n) ~ 1, mu2 ~ 0 + offset(zero),
     theta1 ~ 0 + offset(zero))
}

bcm_exams2_formula <- function() {
  bf(k | trials(n) ~ 1 + (1 | id), mu2 ~ 0 + offset(zero),
     theta1 ~ 0 + offset(zero))
}

test_that("Exams_1 is a two-component mixture with one rate pinned", {
  d <- bcm_exams_data()
  fit <- suppressWarnings(
    frm(bcm_exams1_formula(), family = mixture(binomial, binomial),
        data = d))
  phi <- plogis(unname(fixef(fit)$mu1))
  # the second group knows something, and the first is at chance
  expect_gt(phi, 0.5)
  expect_lt(phi, 1)
  # the pinned pieces have no coefficients at all
  expect_length(fixef(fit)$mu2, 0L)
  expect_length(fixef(fit)$theta1, 0L)
  # and the latent probabilities split the class exactly where the book
  # splits it: the first five people are the guessers
  pr <- mixture_probs(fit)
  expect_true(all(pr[1:5, 2] > 0.5))
  expect_true(all(pr[6:15, 1] > 0.5))
})

test_that("Exams_1 matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_exams_data()
  fit <- suppressWarnings(
    frm(bcm_exams1_formula(), family = mixture(binomial, binomial),
        data = d))
  stan_lp_check(
    bcm_exams1_code(),
    data = list(p = nrow(d), k = as.integer(d$k), n = 40L),
    fit = fit,
    pars = function(f) list(phi = plogis(unname(fixef(f)$mu1))),
    const = 0)
})

test_that("Exams_2 gives the second group individual rates", {
  d <- bcm_exams_data()
  fit <- suppressWarnings(
    frm(bcm_exams2_formula(), family = mixture(binomial, binomial),
        data = d))
  expect_equal(nrow(ranef(fit)[[1L]]), 15L)
  # a random effect with ONE observation per level is an overdispersion
  # term, so the group standard deviation is the whole of what
  # individual differences buy here
  expect_gt(sqrt(unname(VarCorr(fit)[[1L]])[1, 1]), 0)
})

test_that("Exams_2 matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_exams_data()
  fit <- suppressWarnings(
    frm(bcm_exams2_formula(), family = mixture(binomial, binomial),
        data = d))
  stan_lp_check(
    bcm_exams2_code(),
    data = list(p = nrow(d), k = as.integer(d$k), n = 40L),
    fit = fit,
    pars = function(f) {
      b <- unname(fixef(f)$mu1)
      list(mu = b, sigma = sqrt(unname(VarCorr(f)[[1L]])[1, 1]),
           philogit = b + ranef(f)[[1L]][, 1L])
    },
    # The identity compares the JOINT density at the conditional modes,
    # which is exact. What is approximate here is the MARGINAL frmtmb
    # maximizes: the integrand is a mixture of binomials rather than
    # anything Gaussian, so the Laplace approximation to the integral
    # over philogit carries an error the identity cannot see. See the
    # vignette's Laplace caveat, and check_laplace() in frmtmb.sample.
    inner = "philogit",
    const = 0)
})

# STARTING VALUES. mixture() initializes each component's mu from a
# quantile of the response, which for a binomial or beta-binomial is a
# COUNT and falls outside the logit link's range. Both components then
# start at the link's origin, which is the exchangeable point where the
# two are the same distribution, and a beta-binomial mixture stays
# there: the fit comes back with the components identical and every
# latent probability at one half. So every mixture below whose
# components are both free is started from two separated rates, and the
# separation is the only thing the start asserts.
bcm_mix_start <- function(mu, phi = NULL, theta = NULL) {
  out <- list(beta = stats::qlogis(mu))
  bd <- c(if (!is.null(phi)) log(phi), if (!is.null(theta)) theta)
  if (length(bd)) out$betad <- bd
  out
}

# ---------------------------------------------------------------------
# Malingering_1: two groups, both rates free.
# ---------------------------------------------------------------------

bcm_malingering_data <- function() {
  data.frame(k = c(45, 45, 44, 45, 44, 45, 45, 45, 45, 45, 30,
                   20, 6, 44, 44, 27, 25, 17, 14, 27, 35, 30),
             n = 45L, zero = 0)
}

bcm_malingering1_code <- function() {
  paste(
    "data {",
    "  int<lower=1> p; array[p] int<lower=0> k; int<lower=1> n;",
    "}",
    "parameters {",
    "  vector<lower=0, upper=1>[2] psi;",
    "}",
    "model {",
    "  for (i in 1 : p)",
    "    target += log_mix(0.5, binomial_lpmf(k[i] | n, psi[1]),",
    "                           binomial_lpmf(k[i] | n, psi[2]));",
    "}",
    sep = "\n")
}

test_that("Malingering_1 separates two response rates", {
  d <- bcm_malingering_data()
  fit <- suppressWarnings(
    frm(bf(k | trials(n) ~ 1, theta1 ~ 0 + offset(zero)),
        family = mixture(binomial, binomial), data = d,
        start = bcm_mix_start(c(0.55, 0.95))))
  psi <- sort(plogis(c(unname(fixef(fit)$mu1), unname(fixef(fit)$mu2))))
  # the bona fide group answers nearly everything, the malingerers do
  # not; the book's order restriction psi[2] < psi[1] is a PRIOR, and
  # under maximum likelihood the two components are exchangeable, so
  # which label lands on which rate is not identified. sort() is the
  # honest reading of the estimate.
  expect_lt(psi[1], 0.8)
  expect_gt(psi[2], 0.9)
})

test_that("Malingering_1 matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_malingering_data()
  fit <- suppressWarnings(
    frm(bf(k | trials(n) ~ 1, theta1 ~ 0 + offset(zero)),
        family = mixture(binomial, binomial), data = d,
        start = bcm_mix_start(c(0.55, 0.95))))
  stan_lp_check(
    bcm_malingering1_code(),
    data = list(p = nrow(d), k = as.integer(d$k), n = 45L),
    fit = fit,
    pars = function(f) {
      list(psi = plogis(c(unname(fixef(f)$mu1), unname(fixef(f)$mu2))))
    },
    const = 0)
})

# ---------------------------------------------------------------------
# Malingering_2 and Cheating: the same two-group model with individual
# differences INSIDE each group.
#
# The book gives every person a rate of their own with a beta prior
# whose mean and precision belong to their group. That rate appears
# nowhere else, and a binomial rate with a beta prior integrates in
# closed form to a BETA-BINOMIAL. So the model with the per-person rates
# summed out is exactly
#
#   mixture(beta_binomial, beta_binomial)
#
# and frmtmb's beta_binomial is already in the mean-precision
# parameterization the book uses: alpha = mu phi, beta = (1 - mu) phi.
# The integral is not a convenience, it is what makes the model
# estimable at all: left free, one rate per person with one observation
# each is a saturated fit.
# ---------------------------------------------------------------------

bcm_betamix_code <- function() {
  paste(
    "data {",
    "  int<lower=1> p; array[p] int<lower=0> k; int<lower=1> n;",
    "}",
    "parameters {",
    "  real<lower=0, upper=1> w;",
    "  real<lower=0, upper=1> mu1; real<lower=0> phi1;",
    "  real<lower=0, upper=1> mu2; real<lower=0> phi2;",
    "}",
    "model {",
    "  for (i in 1 : p)",
    "    target += log_mix(",
    "      w,",
    "      beta_binomial_lpmf(k[i] | n, mu1 * phi1, (1 - mu1) * phi1),",
    "      beta_binomial_lpmf(k[i] | n, mu2 * phi2, (1 - mu2) * phi2));",
    "}",
    sep = "\n")
}

bcm_betamix_pars <- function(f) {
  g <- function(nm) unname(fixef(f)[[nm]])
  list(w = plogis(g("theta1")),
       mu1 = plogis(g("mu1")), phi1 = exp(g("phi1")),
       mu2 = plogis(g("mu2")), phi2 = exp(g("phi2")))
}

test_that("Malingering_2 is a mixture of beta-binomials", {
  d <- bcm_malingering_data()
  fit <- suppressWarnings(
    frm(bf(k | trials(n) ~ 1),
        family = mixture(beta_binomial, beta_binomial), data = d,
        start = bcm_mix_start(c(0.55, 0.95), phi = c(15, 15), theta = 0)))
  mu <- sort(plogis(c(unname(fixef(fit)$mu1), unname(fixef(fit)$mu2))))
  # the malingerers answer about half, the bona fide group nearly all
  expect_lt(mu[1], 0.7)
  expect_gt(mu[2], 0.9)
  # The first group is genuinely overdispersed and the second is NOT.
  # phi2's profile is monotone up to the binomial limit, so its maximum
  # likelihood estimate is INFINITY, and the fit stops at 7.19e8 only
  # because the surface is flat to machine precision long before that
  # matters. Asserting that the two precisions are "finite" would say
  # nothing at all: fixef() returns log(phi), so such an assertion
  # passes at log(phi2) = 20.4 and cannot fail.
  expect_lt(unname(fixef(fit)$phi1), log(100))
  expect_gt(unname(fixef(fit)$phi2), log(1e6))
  # and this is what an infinite phi2 MEANS, from inside the grammar:
  # the same mixture with component 2 AT the binomial limit reaches the
  # same optimum with one parameter fewer
  lim <- suppressWarnings(
    frm(bf(k | trials(n) ~ 1),
        family = mixture(beta_binomial, binomial), data = d,
        start = bcm_mix_start(c(0.55, 0.95), phi = 15, theta = 0)))
  expect_equal(as.numeric(logLik(lim)), as.numeric(logLik(fit)),
               tolerance = 1e-4)
  expect_equal(sort(plogis(c(unname(fixef(lim)$mu1),
                             unname(fixef(lim)$mu2)))),
               mu, tolerance = 1e-4)
})

test_that("Malingering_2 matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_malingering_data()
  fit <- suppressWarnings(
    frm(bf(k | trials(n) ~ 1),
        family = mixture(beta_binomial, beta_binomial), data = d,
        start = bcm_mix_start(c(0.55, 0.95), phi = c(15, 15), theta = 0)))
  stan_lp_check(
    bcm_betamix_code(),
    data = list(p = nrow(d), k = as.integer(d$k), n = 45L),
    fit = fit, pars = bcm_betamix_pars, const = 0)
})

bcm_cheating_data <- function() {
  k <- c(26, 23, 39, 34, 28, 23, 36, 32, 29, 31, 40, 23, 40, 39, 25, 40,
         34, 29, 38, 37, 33, 30, 37, 31, 40, 34, 38, 23, 39, 19, 35, 26,
         38, 18, 31, 37, 40, 34, 33, 40, 36, 35, 21, 37, 23, 40, 39, 37,
         36, 34, 31, 32, 34, 30, 40, 24, 33, 40, 40, 23, 39, 22, 23, 32,
         37, 40, 25, 34, 40, 27, 35, 32, 40, 36, 40, 29, 30, 35, 39, 28,
         27, 40, 32, 40, 25, 30, 26, 37, 17, 38, 30, 33, 36, 33, 27, 38,
         34, 40, 20, 40, 25, 38, 20, 37, 32, 40, 21, 40, 34, 40, 31, 37,
         27, 34, 36, 36, 21, 26)
  truth <- c(1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1,
             0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1,
             1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0,
             1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
             0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1,
             0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
  data.frame(k = k, n = 40L, truth = truth)
}

test_that("Cheating classifies better than chance", {
  d <- bcm_cheating_data()
  fit <- suppressWarnings(
    frm(bf(k | trials(n) ~ 1),
        family = mixture(beta_binomial, beta_binomial), data = d,
        start = bcm_mix_start(c(0.72, 0.92), phi = c(30, 30), theta = 0)))
  pr <- mixture_probs(fit)
  # the higher-scoring component is the cheating one
  hi <- which.max(c(unname(fixef(fit)$mu1), unname(fixef(fit)$mu2)))
  called <- as.integer(pr[, hi] > 0.5)
  acc <- mean(called == d$truth)
  # the book reports about 60 percent correct; anything at chance would
  # mean the mixture found nothing
  expect_gt(acc, 0.55)
})

# ---------------------------------------------------------------------
# TwentyQuestions: ten people, twenty questions, and one rate for each.
# An answer is correct with probability p_i q_j, a PRODUCT of two rates,
# so this is a nonlinear formula with an identity link rather than a
# two-way additive model on any link scale.
# ---------------------------------------------------------------------

bcm_twenty_matrix <- function() {
  matrix(c(1,1,1,1,0,0,1,1,0,1,0,0,1,0,0,1,0,1,0,0,
           0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,
           0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,
           1,0,1,1,0,1,1,1,0,1,0,0,1,0,0,0,0,1,0,0,
           1,1,0,1,0,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,
           0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,1,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,1,0,1,
           1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0),
         nrow = 10, ncol = 20, byrow = TRUE)
}

bcm_twenty_data <- function() {
  k <- bcm_twenty_matrix()
  data.frame(k = as.vector(k),
             person = factor(rep(seq_len(nrow(k)), times = ncol(k))),
             question = factor(rep(seq_len(ncol(k)), each = nrow(k))))
}

bcm_twenty_code <- function() {
  paste(
    "data {",
    "  int<lower=1> np; int<lower=1> nq;",
    "  array[np, nq] int<lower=0, upper=1> k;",
    "}",
    "parameters {",
    "  vector<lower=0, upper=1>[np] p;",
    "  vector<lower=0, upper=1>[nq] q;",
    "}",
    "model {",
    "  for (i in 1 : np)",
    "    for (j in 1 : nq)",
    "      target += bernoulli_lpmf(k[i, j] | p[i] * q[j]);",
    "}",
    sep = "\n")
}

bcm_twenty_formula <- function() {
  bf(k ~ plogis(lp) * plogis(lq), lp ~ 0 + person, lq ~ 0 + question,
     nl = TRUE)
}

test_that("TwentyQuestions is a product of two rates", {
  d <- bcm_twenty_data()
  fit <- frm(bcm_twenty_formula(),
             family = bernoulli(link = "identity"), data = d)
  p <- plogis(frm_b(fit, "lp"))
  q <- plogis(frm_b(fit, "lq"))
  expect_length(p, 10L)
  expect_length(q, 20L)
  # person 8 got nothing right, so their rate is driven to the boundary,
  # and question 8's rate is the highest of the easy questions
  expect_lt(p[8], min(p[-8]))
  expect_true(all(q >= 0 & q <= 1))
  # the fitted probability of every cell is the product
  expect_equal(unname(fitted(fit)),
               p[as.integer(d$person)] * q[as.integer(d$question)],
               tolerance = 1e-6)
})

test_that("TwentyQuestions matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_twenty_data()
  fit <- frm(bcm_twenty_formula(),
             family = bernoulli(link = "identity"), data = d)
  stan_lp_check(
    bcm_twenty_code(),
    data = list(np = 10L, nq = 20L, k = bcm_twenty_matrix()),
    fit = fit,
    pars = function(f) list(p = plogis(frm_b(f, "lp")),
                            q = plogis(frm_b(f, "lq"))),
    const = 0,
    # A rate driven to a boundary has a gradient the optimizer stops at
    # rather than zeroes, so the value alone is asserted here and the
    # test above checks the map by rebuilding every fitted probability.
    grad = "none")
})

# ---------------------------------------------------------------------
# TwoCountryQuiz: the model the Stan port does not implement.
#
# See inst/bcm/marginal.R for why the double mixture does not factorize
# per row and how the sum is arranged so that it costs 2^questions
# rather than 2^(people + questions).
#
# The reference port carries no data, because it carries no model, so
# the data below has the structure the chapter describes rather than the
# chapter's own numbers: two groups of people, two groups of questions,
# and accuracy that depends on whether they match.
# ---------------------------------------------------------------------

bcm_tcq_matrix <- function() {
  # two blocks, with one miss inside a block and one lucky hit across
  # them, so that neither rate is driven to a boundary
  matrix(c(1, 1, 1, 1, 0, 1, 0, 0,
           1, 1, 1, 1, 0, 0, 0, 0,
           0, 1, 1, 1, 0, 0, 0, 0,
           1, 1, 1, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 1, 1, 1, 1,
           0, 0, 0, 0, 1, 1, 1, 1,
           0, 0, 0, 0, 1, 1, 1, 1,
           0, 0, 0, 0, 1, 1, 1, 0),
         nrow = 8, ncol = 8, byrow = TRUE)
}

bcm_tcq_data <- function() {
  k <- bcm_tcq_matrix()
  data.frame(k = as.vector(k),
             person = factor(rep(seq_len(nrow(k)), times = ncol(k))),
             question = factor(rep(seq_len(ncol(k)), each = nrow(k))))
}

test_that("TwoCountryQuiz recovers the two blocks", {
  skip_unless_bcm("marginal.R")
  d <- bcm_tcq_data()
  fit <- frm(k ~ 1, family = bcm_two_country(person, question), data = d)
  a <- plogis(unname(fixef(fit)$alpha))
  b <- plogis(unname(fixef(fit)$beta))
  # own-country questions are answered, other-country ones are not
  expect_gt(a, 0.8)
  expect_lt(b, 0.2)
  # Which country is called 1 is NOT identified: every question pattern
  # has a complement under which every person's country flips, and the
  # two carry the same likelihood. So every person's marginal is exactly
  # one half, and that is the right answer rather than a failure.
  pr <- latent_probs(fit)
  expect_equal(dim(pr), c(8L, 2L))
  expect_equal(unname(pr), matrix(0.5, 8L, 2L), tolerance = 1e-8)
  # The quantity that survives the flip is whether two people came from
  # the SAME country, and that recovers the blocks exactly.
  ag <- bcm_two_country_agreement(fit)
  expect_equal(dim(ag), c(8L, 8L))
  same <- outer(rep(1:2, each = 4), rep(1:2, each = 4), `==`)
  expect_gt(min(ag[same]), 0.9)
  expect_lt(max(ag[!same]), 0.1)
})

test_that("TwoCountryQuiz refuses a quiz it cannot enumerate", {
  skip_unless_bcm("marginal.R")
  k <- matrix(rbinom(20 * 20, 1L, 0.5), 20, 20)
  d <- data.frame(k = as.vector(k),
                  person = factor(rep(1:20, times = 20)),
                  question = factor(rep(1:20, each = 20)))
  expect_error(
    frm(k ~ 1, family = bcm_two_country(person, question), data = d),
    "exponential in the number of")
})
