# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 19: the
# Balloon Analogue Risk Task.
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/CaseStudies/BART.
#
# The model looks like a cognitive process model and IS a Bernoulli
# logistic regression. Each pump k of a balloon is a decision, and the
# probability of stopping there is
#
#   1 - inv_logit(-beta (k - omega)) = inv_logit(beta k - beta omega),
#
# so the linear predictor is beta * k with an intercept of -beta omega.
# Fitting `stop ~ k` recovers both, and the book's two parameters come
# back by arithmetic: omega = -intercept / slope is the optimal number
# of pumps, and gamma+ = -omega log(1 - p) is the risk-taking parameter
# behind it. Nothing here needs a custom family.
#
# The data is George's, from GeorgeSober.txt, GeorgeTipsy.txt and
# GeorgeDrunk.txt: 90 balloons per condition, reduced to how many pumps
# he made and whether he cashed out or the balloon burst.

bcm_bart_trials <- function() {
  list(
    sober = list(
      pumps = c(2, 1, 3, 3, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 5, 2, 4, 4,
                4, 4, 2, 1, 2, 3, 1, 1, 1, 1, 1, 5, 5, 4, 1, 3, 1, 3, 3,
                4, 4, 6, 2, 2, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 7, 4, 4, 3,
                1, 1, 1, 2, 3, 2, 3, 2, 3, 2, 4, 2, 5, 1, 1, 3, 2, 5, 1,
                3, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1),
      cash = c(0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,
               1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1,
               1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
               1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0,
               1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1)),
    tipsy = list(
      pumps = c(4, 4, 2, 3, 5, 5, 5, 4, 1, 3, 3, 4, 1, 1, 4, 1, 1, 1, 5,
                3, 5, 5, 2, 2, 6, 7, 5, 1, 5, 4, 6, 5, 4, 5, 5, 5, 5, 5,
                5, 2, 5, 5, 5, 5, 7, 1, 2, 1, 5, 5, 5, 4, 3, 7, 7, 7, 8,
                7, 8, 7, 1, 4, 4, 1, 4, 2, 4, 4, 5, 4, 4, 1, 4, 4, 1, 4,
                1, 1, 2, 4, 4, 4, 2, 2, 4, 3, 4, 3, 4, 3),
      cash = c(1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1,
               1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
               1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1,
               1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1,
               0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0)),
    drunk = list(
      pumps = c(1, 7, 7, 6, 5, 4, 3, 6, 4, 4, 2, 1, 4, 4, 2, 3, 3, 3, 3,
                6, 5, 5, 5, 2, 5, 2, 4, 5, 5, 2, 9, 1, 6, 6, 6, 6, 4, 6,
                5, 5, 2, 6, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 3, 5, 1, 2, 2, 2, 3, 3, 3, 2, 2, 2, 2, 2, 3, 3,
                3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1),
      cash = c(1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1,
               1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1,
               0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
               1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1,
               1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1)))
}

# One row per DECISION. A balloon pumped npumps times and then cashed
# offers npumps + 1 decisions, the last of which was to stop; a balloon
# that burst offers npumps, all of them decisions to keep going.
bcm_bart_long <- function(conds = "sober") {
  tr <- bcm_bart_trials()
  out <- list()
  for (cn in conds) {
    z <- tr[[cn]]
    opts <- z$pumps + z$cash
    k <- unlist(lapply(opts, seq_len))
    stop_ <- unlist(Map(function(o, c_) c(rep(0L, o - c_), rep(1L, c_)),
                        opts, z$cash))
    out[[cn]] <- data.frame(stop = stop_, k = k, cond = cn,
                            trial = rep(seq_along(opts), opts))
  }
  d <- do.call(rbind, out)
  d$cond <- factor(d$cond, levels = conds)
  rownames(d) <- NULL
  d
}

bcm_bart_code <- function() {
  paste(
    "data {",
    "  int<lower=1> N;",
    "  array[N] int<lower=0, upper=1> d;",
    "  vector[N] k;",
    "  real<lower=0, upper=1> p;",
    "}",
    "parameters {",
    "  real<lower=0> gplus;",
    "  real<lower=0> beta;",
    "}",
    "model {",
    "  real omega = -gplus / log1m(p);",
    "  for (i in 1 : N)",
    "    target += bernoulli_lpmf(d[i] |",
    "                             1 - inv_logit(-beta * (k[i] - omega)));",
    "}",
    sep = "\n")
}

# The book's two parameters, from the two regression coefficients.
bcm_bart_pars <- function(fit, p = 0.15, dpar = "mu",
                          icept = "(Intercept)", slope = "k") {
  b <- fixef(fit)[[dpar]]
  beta <- unname(b[slope])
  omega <- -unname(b[icept]) / beta
  list(beta = beta, omega = omega, gplus = -omega * log1p(-p))
}

# ---------------------------------------------------------------------
# BART_1: one condition.
# ---------------------------------------------------------------------

test_that("BART_1 is a Bernoulli logit regression on the pump index", {
  d <- bcm_bart_long("sober")
  fit <- frm(stop ~ k, family = bernoulli(), data = d)
  pp <- bcm_bart_pars(fit)
  # a positive slope is a person who becomes more likely to stop the
  # further they pump, which is what the model says everybody does
  expect_gt(pp$beta, 0)
  # the optimal number of pumps is in the range the balloons allow
  expect_gt(pp$omega, 0)
  expect_lt(pp$omega, 10)
  # and the risk parameter is positive, which is the book's gamma+
  expect_gt(pp$gplus, 0)
})

test_that("BART_1 matches its Stan program", {
  skip_unless_stan_identity()
  d <- bcm_bart_long("sober")
  fit <- frm(stop ~ k, family = bernoulli(), data = d)
  stan_lp_check(
    bcm_bart_code(),
    data = list(N = nrow(d), d = as.integer(d$stop), k = as.numeric(d$k),
                p = 0.15),
    fit = fit,
    pars = function(f) {
      pp <- bcm_bart_pars(f)
      list(gplus = pp$gplus, beta = pp$beta)
    },
    # the original's <lower=0, upper=10> on both parameters is a uniform
    # prior it does not add to the target, and the estimates are inside
    # it, so the two programs hold the same function
    const = 0)
})

# ---------------------------------------------------------------------
# BART_2: three conditions.
#
# The original draws each condition's pair from a normal truncated at
# zero. With THREE levels that group distribution is a prior rather than
# something the data identifies, so the fit here is the unpooled model,
# which is the likelihood the original's prior multiplies. The book's
# own conclusion, that risk taking rises with intoxication, is a
# statement about the three estimates and survives unchanged.
# ---------------------------------------------------------------------

test_that("BART_2 gives every condition its own risk parameter", {
  d <- bcm_bart_long(c("sober", "tipsy", "drunk"))
  fit <- frm(stop ~ 0 + cond + cond:k, family = bernoulli(), data = d)
  g <- vapply(c("sober", "tipsy", "drunk"), function(cn) {
    bcm_bart_pars(fit, icept = paste0("cond", cn),
                  slope = paste0("cond", cn, ":k"))$gplus
  }, 0)
  expect_length(g, 3L)
  expect_true(all(g > 0))
  # the book's finding: George takes more risk the drunker he is
  expect_lt(g[["sober"]], g[["tipsy"]])
  expect_lt(g[["tipsy"]], g[["drunk"]])
})
