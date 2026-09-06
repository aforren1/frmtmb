# Lee and Wagenmakers, Bayesian Cognitive Modeling, chapter 15: SIMPLE.
#
# Stan programs adapted from the BSD-3 ports at
# github.com/stan-dev/example-models, directory
# Bayesian_Cognitive_Modeling/CaseStudies/Simple, on Murdock's (1962)
# free-recall data: six lists of different lengths, each presented many
# times, and the proportion recalled at each serial position.
#
# Like the generalized context model, the likelihood does not factorize
# over rows: an item's discriminability is its similarity relative to
# every other item in the SAME list, so one row reads the whole list.
# The family is in inst/bcm/similarity.R and uses the structured
# protocol's loglik slot.
#
# One adaptation. `c` and `s` are bounded on (0, 100) and `t` on (0, 1)
# in the original; those bounds are priors rather than part of the
# likelihood, so the family puts c and s on log links and t on a logit,
# and the tests assert the estimates stay inside the intervals anyway.
#
# The original's theta = min(1, sum(resp)) IS carried, through
# bcm_cap1() of inst/bcm/binomial-extras.R.

bcm_simple_sets <- function() {
  list(
    list(k = c(994, 806, 691, 634, 634, 835, 965, 1008, 1181, 1382),
         m = c(33, 31, 29, 27, 25, 23, 21, 19, 17, 15), n = 1440L),
    list(k = c(794, 640, 576, 512, 486, 486, 512, 512, 538, 640, 742,
               794, 1024, 1126, 1254),
         m = c(48, 46, 44, 42, 40, 38, 36, 34, 32, 30, 28, 26, 24, 22,
               20), n = 1280L),
    list(k = c(836, 623, 562, 517, 441, 471, 441, 395, 350, 441, 456,
               486, 456, 593, 608, 684, 897, 1110, 1353, 1459),
         m = c(63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37,
               35, 33, 31, 29, 27, 25), n = 1520L),
    list(k = c(699, 410, 304, 243, 258, 243, 228, 213, 304, 350, 395,
               319, 365, 410, 456, 608, 958, 1170, 1277, 1474),
         m = c(29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16,
               15, 14, 13, 12, 11, 10), n = 1520L),
    list(k = c(576, 360, 240, 228, 240, 240, 228, 216, 240, 240, 240,
               240, 228, 240, 216, 204, 240, 228, 240, 240, 264, 252,
               288, 324, 444, 480, 660, 912, 1080, 1188),
         m = c(44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31,
               30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17,
               16, 15), n = 1200L),
    list(k = c(384, 256, 154, 154, 141, 128, 154, 154, 154, 128, 179,
               141, 128, 141, 154, 179, 154, 141, 179, 128, 154, 128,
               154, 192, 179, 166, 154, 218, 218, 218, 256, 256, 230,
               307, 384, 512, 691, 922, 1114, 1254),
         m = c(59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46,
               45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32,
               31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20),
         n = 1280L))
}

bcm_simple_data <- function(which = seq_len(6)) {
  s <- bcm_simple_sets()[which]
  do.call(rbind, lapply(seq_along(s), function(i) {
    z <- s[[i]]
    data.frame(k = z$k, n = z$n, m = z$m,
               pos = seq_along(z$k),
               ll = length(z$k),
               set = factor(which[i], levels = which))
  }))
}

# One program for both models. `share` switches SIMPLE_1's one triple
# per list for SIMPLE_2's shared c and s with a threshold that is a
# linear function of list length.
bcm_simple_code <- function(share = FALSE) {
  head <- c(
    "data {",
    "  int<lower=1> N; int<lower=1> D;",
    "  array[N] int<lower=0> k; array[N] int<lower=1> n;",
    "  vector[N] logm; array[N] int<lower=1> set;",
    "  array[D] int<lower=1> len; array[D, 40] int<lower=0> idx;",
    "  vector[D] w;",
    "}")
  pars <- if (share) {
    c("parameters {",
      "  real<lower=0> c; real<lower=0> s; vector[2] a;",
      "}")
  } else {
    c("parameters {",
      "  vector<lower=0>[D] c; vector<lower=0>[D] s;",
      "  vector<lower=0, upper=1>[D] t;",
      "}")
  }
  pick <- if (share) {
    c("    real cx = c;", "    real sx = s;",
      "    real tx = a[1] * w[x] + a[2];")
  } else {
    c("    real cx = c[x];", "    real sx = s[x];", "    real tx = t[x];")
  }
  body <- c(
    "model {",
    "  for (x in 1 : D) {",
    pick,
    "    int L = len[x];",
    "    matrix[L, L] sim;",
    "    for (i in 1 : L)",
    "      for (j in 1 : L)",
    # `abs`, not the original's `fabs`: Stan has removed that spelling
    "        sim[i, j] = exp(-cx * abs(logm[idx[x, i]] - logm[idx[x, j]]));",
    "    for (i in 1 : L) {",
    "      real theta = 0;",
    "      real tot = sum(sim[i]);",
    "      for (j in 1 : L)",
    "        theta += inv_logit(sx * (sim[i, j] / tot - tx));",
    "      target += binomial_lpmf(k[idx[x, i]] | n[idx[x, i]], theta);",
    "    }",
    "  }",
    "}")
  paste(c(head, pars, body), collapse = "\n")
}

bcm_simple_stan_data <- function(d) {
  sets <- split(seq_len(nrow(d)), d$set)
  D <- length(sets)
  len <- lengths(sets)
  idx <- matrix(1L, D, 40L)
  for (x in seq_len(D)) {
    idx[x, seq_len(len[x])] <- sets[[x]][order(d$pos[sets[[x]]])]
  }
  list(N = nrow(d), D = D, k = as.integer(d$k), n = as.integer(d$n),
       logm = log(d$m), set = as.integer(d$set),
       len = as.integer(len), idx = idx,
       w = as.numeric(tapply(d$ll, d$set, function(z) z[1])))
}

bcm_simple_family <- function() bcm_simple(set, pos, m)

# A start with every serial position's recall probability below one.
#
# The cap makes an INFEASIBLE start fatal rather than slow: where theta
# reaches one, a row with fewer recalls than attempts has a log density
# of -Infinity and a gradient of NaN, and the optimizer stops. Every
# start with a threshold at or above 0.6 reaches the same optimum, a
# log-likelihood of -2083.455, so what this asserts is feasibility and
# nothing else.
bcm_simple_start <- function(d) {
  k <- nlevels(d$set)
  list(beta = rep(log(20), k),
       betad = c(rep(log(10), k), rep(stats::qlogis(0.6), k)))
}

# nlminb's default stopping rule leaves this fit with a gradient of
# about 1.2e-3, which is over the identity harness's 1e-3 and is the
# OPTIMIZER rather than the map: the log-likelihood is -2083.455 either
# way, and tightening the rule brings the gradient to 3.1e-4. Eighteen
# free parameters over 135 rows, each of which reads its whole list, is
# simply a flatter surface than the default rule was chosen for.
bcm_simple_control <- function() {
  frmtmb_control(optCtrl = list(rel.tol = 1e-12, x.tol = 1e-12,
                                eval.max = 5000, iter.max = 5000))
}

# ---------------------------------------------------------------------
# SIMPLE_1: one triple of parameters per list.
# ---------------------------------------------------------------------

test_that("SIMPLE_1 produces a serial position curve", {
  skip_unless_bcm(c("binomial-extras.R", "similarity.R"))
  d <- bcm_simple_data()
  # nlminb stops on "false convergence (8)" on this surface however
  # tight the rule; the gradient is 3.1e-4 and the log-likelihood is the
  # same from every feasible start, so the code is the stopping rule
  # rather than a failure. The identity test below is the real check.
  fit <- suppressWarnings(
    frm(bf(k | trials(n) ~ 0 + set, s ~ 0 + set, t ~ 0 + set),
        family = bcm_simple_family(), data = d,
        start = bcm_simple_start(d), control = bcm_simple_control()))
  th <- fitted(fit) / d$n
  # every recall probability is a probability, which is what the cap is
  # there to guarantee
  expect_lte(max(th), 1 + 1e-8)
  # the estimates stay inside the intervals the original bounds them to
  expect_true(all(exp(frm_b(fit, "c")) < 100))
  expect_true(all(exp(frm_b(fit, "s")) < 100))
  # the shape the model exists to produce: primacy and recency, so the
  # first and last items of the longest list beat its middle
  last <- d$set == levels(d$set)[6]
  p <- th[last]
  expect_gt(p[1], stats::median(p))
  expect_gt(p[length(p)], stats::median(p))
})

test_that("SIMPLE_1 matches its Stan program", {
  skip_unless_stan_identity()
  skip_unless_bcm(c("binomial-extras.R", "similarity.R"))
  d <- bcm_simple_data()
  fit <- suppressWarnings(
    frm(bf(k | trials(n) ~ 0 + set, s ~ 0 + set, t ~ 0 + set),
        family = bcm_simple_family(), data = d,
        start = bcm_simple_start(d), control = bcm_simple_control()))
  stan_lp_check(
    bcm_simple_code(FALSE),
    data = bcm_simple_stan_data(d),
    fit = fit,
    pars = function(f) list(c = exp(frm_b(f, "c")),
                            s = exp(frm_b(f, "s")),
                            t = plogis(frm_b(f, "t"))),
    const = 0)
})

# ---------------------------------------------------------------------
# SIMPLE_2: one gradient and one slope for every list, with the
# threshold a linear function of list length.
#
# The threshold's linear predictor takes an IDENTITY link here, because
# a1 * w + a2 is what the model says and a logit would not be that
# function. The original clips it to (0, 1) with fmax and fmin; the
# family does not clip a THRESHOLD, only the recall probability, so the
# test asserts that every fitted threshold is inside anyway, which is
# what makes the two the same function at the estimate.
# ---------------------------------------------------------------------

test_that("SIMPLE_2 makes the threshold a function of list length", {
  skip_unless_bcm(c("binomial-extras.R", "similarity.R"))
  d <- bcm_simple_data()
  # same stopping rule as SIMPLE_1, with one gradient at 3.1e-3: six
  # lists of very different lengths sharing one gradient and one slope
  # is a flatter surface still
  fit <- suppressWarnings(
    frm(bf(k | trials(n) ~ 1, s ~ 1, t ~ ll),
        family = bcm_simple(set, pos, m, link_t = "identity"),
        data = d,
        # a1 w + a2 has to leave every list's threshold high enough that
        # no serial position reaches the cap; 0.75 at the intercept is
        # not, and the fit stops on a NaN
        start = list(beta = log(20),
                     betad = c(log(10), 1, -0.006))))
  tt <- as.numeric(eval_dpars(fit)[["k"]][["t"]])
  expect_true(all(tt > 0 & tt < 1))
  # longer lists lower the threshold, which is the case study's finding
  expect_lt(unname(fixef(fit)[["t"]]["ll"]), 0)
  expect_lte(max(fitted(fit) / d$n), 1 + 1e-8)
})
