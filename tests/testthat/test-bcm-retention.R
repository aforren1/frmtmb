# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 10: memory
# retention.
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/CaseStudies/MemoryRetention. Sampling
# statements become `target += ..._lpmf(...)`, and the padded matrix
# with its excluded fourth subject and tenth lag becomes one long data
# frame with those cells simply absent.
#
# The retention function is theta = min(1, exp(-alpha t) + beta), a
# NONLINEAR predictor with an identity link, so the frmtmb spelling is
# an nl formula. alpha and beta live on (0, 1), which the formula
# imposes with plogis() rather than a link, because a nonlinear
# parameter has no link of its own.
#
# The `min(1, .)` cap is carried, through bcm_cap1() of
# inst/bcm/binomial-extras.R. It is not decoration: subject 1 recalled
# every item at the shortest lag, so without a cap the fit drives
# exp(-alpha t) + beta past one and the likelihood of that cell grows
# without bound. Retention_1 does not reach the cap and Retention_2 and
# Retention_3 do, which is why all three are fitted the same way.

bcm_retention_data <- function() {
  t <- c(1, 2, 4, 7, 12, 21, 35, 59, 99)
  k <- c(18, 18, 16, 13, 9, 6, 4, 4, 4,
         17, 13, 9, 6, 4, 4, 4, 4, 4,
         14, 10, 6, 4, 4, 4, 4, 4, 4)
  data.frame(k = k, n = 18L, t = rep(t, 3),
             id = factor(rep(1:3, each = length(t))))
}

# One program for all three, with the group structure switched by data:
#   np = 1  one alpha and one beta for everybody (Retention_1)
#   np = 2  one pair per subject, unpooled          (Retention_2)
#   np = 3  one pair per subject from a group distribution (Retention_3)
#
# Retention_3's group distribution is the ONE adaptation of substance in
# this file. The original draws alpha_i and beta_i from normals
# truncated to (0, 1); frmtmb's random effects are normal on the linear
# predictor, which here is the logit, so the group distribution is
# logit-normal. That is a different prior for the same model and it is
# stated as such: a truncated-normal random effect is not in this
# grammar, and nothing in the port pretends otherwise.
bcm_retention_code <- function(np) {
  head <- c(
    "data {",
    "  int<lower=1> N; int<lower=1> S;",
    "  array[N] int<lower=0> k; int<lower=1> n;",
    "  vector[N] t; array[N] int<lower=1> subj;",
    "}")
  pars <- switch(
    as.character(np),
    "1" = c("parameters {",
            "  real<lower=0, upper=1> alpha;",
            "  real<lower=0, upper=1> beta;",
            "}"),
    "2" = c("parameters {",
            "  vector<lower=0, upper=1>[S] alpha;",
            "  vector<lower=0, upper=1>[S] beta;",
            "}"),
    "3" = c("parameters {",
            "  real alphamu; real betamu;",
            "  real<lower=0> alphasigma; real<lower=0> betasigma;",
            "  vector[S] alphalogit; vector[S] betalogit;",
            "}"))
  body <- switch(
    as.character(np),
    "1" = c("  vector[N] theta;",
            "  for (i in 1 : N)",
            "    theta[i] = fmin(1, exp(-alpha * t[i]) + beta);"),
    "2" = c("  vector[N] theta;",
            "  for (i in 1 : N)",
            "    theta[i] = fmin(1, exp(-alpha[subj[i]] * t[i])",
            "                       + beta[subj[i]]);"),
    "3" = c("  vector[S] alpha = inv_logit(alphalogit);",
            "  vector[S] beta = inv_logit(betalogit);",
            "  vector[N] theta;",
            "  target += normal_lpdf(alphalogit | alphamu, alphasigma);",
            "  target += normal_lpdf(betalogit | betamu, betasigma);",
            "  for (i in 1 : N)",
            "    theta[i] = fmin(1, exp(-alpha[subj[i]] * t[i])",
            "                       + beta[subj[i]]);"))
  paste(c(head, pars, "model {", body,
          "  target += binomial_lpmf(k | n, theta);", "}"),
        collapse = "\n")
}

# Starting values. A nonlinear formula has no design matrix to read a
# scale off, so frmtmb starts every nonlinear parameter at zero, and
# exp(-plogis(0) * 99) + plogis(0) is 0.5 for a cell whose observed rate
# is 4/18: the gradient there is fine, but the hierarchical version
# reaches an NaN before it gets anywhere. These are the book's own
# neighbourhood, one logit each.
bcm_retention_start <- function(np) {
  la <- stats::qlogis(0.35)
  lb <- stats::qlogis(0.22)
  b <- if (np == 2L) c(rep(la, 3), rep(lb, 3)) else c(la, lb)
  list(beta = b)
}

bcm_retention_stan_data <- function(d) {
  list(N = nrow(d), S = nlevels(d$id), k = as.integer(d$k),
       n = 18L, t = as.numeric(d$t), subj = as.integer(d$id))
}

# ---------------------------------------------------------------------
# Retention_1: one retention function for everybody.
# ---------------------------------------------------------------------

bcm_retention1_formula <- function() {
  bf(k | trials(n) ~ bcm_cap1(exp(-plogis(la) * t) + plogis(lb)),
     la ~ 1, lb ~ 1, nl = TRUE)
}

test_that("Retention_1 is a nonlinear binomial with an identity link", {
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_retention_data()
  fit <- frm(bcm_retention1_formula(),
             family = binomial(link = "identity"), data = d,
             start = bcm_retention_start(1L))
  a <- plogis(unname(fixef(fit)$la))
  b <- plogis(unname(fixef(fit)$lb))
  # forgetting, and a floor the curve does not fall through
  expect_gt(a, 0)
  expect_gt(b, 0)
  expect_lt(b, 1)
  # every fitted retention rate is a probability, which is what the cap
  # is there to guarantee
  expect_lte(max(fitted(fit) / d$n), 1 + 1e-8)
})

test_that("Retention_1 matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_retention_data()
  fit <- frm(bcm_retention1_formula(),
             family = binomial(link = "identity"), data = d,
             start = bcm_retention_start(1L))
  stan_lp_check(
    bcm_retention_code(1L),
    data = bcm_retention_stan_data(d),
    fit = fit,
    pars = function(f) list(alpha = plogis(unname(fixef(f)$la)),
                            beta = plogis(unname(fixef(f)$lb))),
    # beta(1, 1) on each rate is the uniform density on the unit
    # interval, whose log is zero
    const = 0)
})

# ---------------------------------------------------------------------
# Retention_2: one retention function per subject, unpooled.
# ---------------------------------------------------------------------

bcm_retention2_formula <- function() {
  bf(k | trials(n) ~ bcm_cap1(exp(-plogis(la) * t) + plogis(lb)),
     la ~ 0 + id, lb ~ 0 + id, nl = TRUE)
}

test_that("Retention_2 gives every subject its own pair", {
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_retention_data()
  fit <- frm(bcm_retention2_formula(),
             family = binomial(link = "identity"), data = d,
             start = bcm_retention_start(2L))
  expect_length(fixef(fit)$la, 3L)
  expect_length(fixef(fit)$lb, 3L)
  # the book's point: the three subjects differ, and the third forgets
  # fastest
  a <- plogis(unname(fixef(fit)$la))
  expect_equal(which.max(a), 3L)
  expect_lte(max(fitted(fit) / d$n), 1 + 1e-8)
})

test_that("Retention_2 matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_retention_data()
  fit <- frm(bcm_retention2_formula(),
             family = binomial(link = "identity"), data = d,
             start = bcm_retention_start(2L))
  stan_lp_check(
    bcm_retention_code(2L),
    data = bcm_retention_stan_data(d),
    fit = fit,
    pars = function(f) list(alpha = plogis(frm_b(f, "la")),
                            beta = plogis(frm_b(f, "lb"))),
    const = 0)
})

# ---------------------------------------------------------------------
# Retention_3: the subjects' pairs come from group distributions.
# ---------------------------------------------------------------------

bcm_retention3_formula <- function() {
  bf(k | trials(n) ~ bcm_cap1(exp(-plogis(la) * t) + plogis(lb)),
     la ~ 1 + (1 | id), lb ~ 1 + (1 | id), nl = TRUE)
}

test_that("Retention_3 gives the subjects a group distribution", {
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_retention_data()
  # nlminb probes past the cap on its way in and reports the NaN it
  # finds there before backing off; the fit converges. bcm_cap1() in
  # inst/bcm/binomial-extras.R explains the edge.
  fit <- suppressWarnings(
    frm(bcm_retention3_formula(), family = binomial(link = "identity"),
        data = d, start = bcm_retention_start(3L)))
  # both group standard deviations are estimated and finite, which is
  # the whole difference from Retention_2
  expect_true(is.finite(frm_sd_term(fit, "la: 1 | id")))
  expect_true(is.finite(frm_sd_term(fit, "lb: 1 | id")))
  expect_equal(length(frm_u_term(fit, "la: 1 | id")), 3L)
  expect_lte(max(fitted(fit) / d$n), 1 + 1e-8)
  # and the fitted curve still falls with the lag, for every subject
  for (s in levels(d$id)) {
    p <- fitted(fit)[d$id == s] / 18
    expect_lte(p[length(p)], p[1] + 1e-8)
  }
})

test_that("Retention_3 matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm("binomial-extras.R")
  d <- bcm_retention_data()
  fit <- suppressWarnings(
    frm(bcm_retention3_formula(), family = binomial(link = "identity"),
        data = d, start = bcm_retention_start(3L)))
  stan_lp_check(
    bcm_retention_code(3L),
    data = bcm_retention_stan_data(d),
    fit = fit,
    pars = function(f) {
      list(alphamu = unname(fixef(f)$la), betamu = unname(fixef(f)$lb),
           alphasigma = frm_sd_term(f, "la: 1 | id"),
           betasigma = frm_sd_term(f, "lb: 1 | id"),
           alphalogit = unname(fixef(f)$la) + frm_u_term(f, "la: 1 | id"),
           betalogit = unname(fixef(f)$lb) + frm_u_term(f, "lb: 1 | id"))
    },
    inner = c("alphalogit", "betalogit"),
    const = 0)
})
