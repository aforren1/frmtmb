# One ordinal threshold by its number. brms 2.23.0's default_prior()
# lists class "Intercept" once and then once per threshold, coef = "1",
# "2", ..., under flexible thresholds, under each level of thres(gr = )
# when they are grouped; it lists none under equidistant thresholds and
# refuses coef there ("The following priors do not correspond to any
# model parameter"), and it accepts prior(coef = "2") as a density on
# Intercept[2] (dev/fixes-prior-brms.R, dev/fixes-log/prior-brms.txt).
# Up to 0.67.0 frmtmb listed no coef rows and refused the prior.
#
# And confint() names the internal threshold parameters by what they
# are, where it named them tau_raw_<k> up to 0.67.0.

otn_data <- function(seed = 62, n = 200) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
  u <- stats::rlogis(n) + 0.8 * d$x
  d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
  d$yg <- ifelse(d$g == "b", pmin(d$y, 4L), d$y)
  d
}

# the log prior a MAP fit carries, at its estimates
otn_log_prior <- function(fit) {
  ent <- frmtmb:::resolve_prior_input(list(frame = fit$frame,
                                           spec = fit$spec),
                                      fit$prior)$entries
  -frmtmb:::neg_log_prior_fn(ent)(fit$estimates)
}

otn_tau <- function(fit) {
  frmtmb:::ord_threshold_values(family(fit), fit$estimates$tau_raw)
}

test_that("default_prior() lists a row per threshold, as brms does", {
  d <- otn_data()
  rows <- function(dp) {
    r <- dp[dp$class == "Intercept", ]
    paste(r$group, r$coef, sep = "|")
  }
  dp <- default_prior(y ~ x, data = d, family = cumulative())
  expect_identical(rows(dp), c("|", "|1", "|2", "|3", "|4"))
  dp <- default_prior(y ~ x, data = d, family = sratio())
  expect_identical(rows(dp), c("|", "|1", "|2", "|3", "|4"))
  # grouped: under each level, as many rows as that level's thresholds
  dp <- default_prior(yg | thres(gr = g) ~ x, data = d,
                      family = cumulative())
  expect_identical(rows(dp), c("|", "a|", "a|1", "a|2", "a|3", "a|4",
                               "b|", "b|1", "b|2", "b|3"))
  # equidistant: the class row alone
  dp <- default_prior(y ~ x, data = d,
                      family = cumulative(threshold = "equidistant"))
  expect_identical(rows(dp), "|")
})

test_that("a prior on one threshold is that threshold's density", {
  d <- otn_data()
  # no predictor, so the thresholds are not centered and the density is
  # on the threshold itself
  pl <- set_prior("normal(0, 2)", class = "Intercept") +
    set_prior("normal(-0.4, 0.3)", class = "Intercept", coef = "2")
  fit <- frm(y ~ 1, family = cumulative(), data = d, prior = pl)
  tau <- otn_tau(fit)
  raw <- fit$estimates$tau_raw
  # the coef row replaces the class density at threshold 2 only; the
  # ordered map's log-Jacobian enters once
  ref <- sum(stats::dnorm(tau[-2L], 0, 2, log = TRUE)) +
    stats::dnorm(tau[2L], -0.4, 0.3, log = TRUE) + sum(raw[-1L])
  expect_lt(abs(otn_log_prior(fit) - ref) / abs(ref), 1e-10)
  # a coef row alone puts a density on that threshold and nothing else
  fit1 <- frm(y ~ 1, family = cumulative(), data = d,
              prior = set_prior("normal(-0.4, 0.3)", class = "Intercept",
                                coef = "2"))
  tau <- otn_tau(fit1)
  ref <- stats::dnorm(tau[2L], -0.4, 0.3, log = TRUE) +
    sum(fit1$estimates$tau_raw[-1L])
  expect_lt(abs(otn_log_prior(fit1) - ref) / abs(ref), 1e-10)
  # the penalty moves the estimate toward the prior location
  f0 <- frm(y ~ 1, family = cumulative(), data = d)
  t0 <- otn_tau(f0)
  tight <- frm(y ~ 1, family = cumulative(), data = d,
               prior = set_prior("normal(-0.1, 0.001)", class = "Intercept",
                                 coef = "2"))
  expect_lt(abs(otn_tau(tight)[2L] + 0.1), abs(t0[2L] + 0.1) / 10)
  # an unordered family holds the thresholds themselves: no Jacobian
  fs <- frm(y ~ 1, family = sratio(), data = d,
            prior = set_prior("normal(0.3, 0.5)", class = "Intercept",
                              coef = "3"))
  ref <- stats::dnorm(otn_tau(fs)[3L], 0.3, 0.5, log = TRUE)
  expect_lt(abs(otn_log_prior(fs) - ref) / abs(ref), 1e-10)
  # one threshold of one level of grouped thresholds
  fg <- frm(yg | thres(gr = g) ~ 1, family = cumulative(), data = d,
            prior = set_prior("normal(0.5, 0.4)", class = "Intercept",
                              group = "b", coef = "3"))
  tau <- otn_tau(fg)
  # level a has 4 thresholds, so b's third is the 7th, and b's slice is
  # internal 5:7
  ref <- stats::dnorm(tau[7L], 0.5, 0.4, log = TRUE) +
    sum(fg$estimates$tau_raw[6:7])
  expect_lt(abs(otn_log_prior(fg) - ref) / abs(ref), 1e-10)
})

test_that("a threshold coef is refused where brms has no such parameter", {
  d <- otn_data()
  one <- function(fam, coef, group = "", f = y ~ x) {
    tryCatch({
      frm(f, family = fam, data = d,
          prior = set_prior("normal(0, 1)", class = "Intercept",
                            coef = coef, group = group))
      NA_character_
    }, error = function(e) conditionMessage(e))
  }
  expect_match(one(cumulative(threshold = "equidistant"), "2"),
               "takes no coef", fixed = TRUE)
  expect_match(one(cumulative(threshold = "sum_to_zero"), "2"),
               "sum_to_zero", fixed = TRUE)
  expect_match(one(cumulative(), "5"), "has 4 thresholds", fixed = TRUE)
  expect_match(one(cumulative(), "x"), "numbers a threshold", fixed = TRUE)
  expect_match(one(cumulative(), "2", f = yg | thres(gr = g) ~ x),
               "group = and coef =", fixed = TRUE)
  expect_match(one(cumulative(), "4", group = "b",
                   f = yg | thres(gr = g) ~ x),
               "has 3 thresholds", fixed = TRUE)
})

test_that("confint() names the thresholds by what they are", {
  d <- otn_data()
  lab <- function(fit) rownames(confint(fit))
  expect_identical(lab(frm(y ~ x, family = sratio(), data = d)),
                   c("x", paste0("Intercept[", 1:4, "]")))
  fc <- frm(y ~ x, family = cumulative(), data = d)
  expect_identical(lab(fc), c("x", "Intercept[1]",
                              "log(Intercept[2] - Intercept[1])",
                              "log(Intercept[3] - Intercept[2])",
                              "log(Intercept[4] - Intercept[3])"))
  expect_identical(rownames(vcov(fc, full = TRUE)), lab(fc))
  # the first threshold IS the internal parameter: the row agrees with
  # fixef() exactly
  expect_identical(unname(confint(fc, parm = "Intercept[1]")[, "est"]),
                   unname(fixef(fc)["Intercept[1]", "Estimate"]))
  expect_identical(confint(fc, parm = "b_Intercept[1]"),
                   confint(fc, parm = "Intercept[1]"))
  fe <- frm(y ~ x, family = cumulative(threshold = "equidistant"), data = d)
  expect_identical(lab(fe), c("x", "Intercept[1]", "log(delta)"))
  ci <- NULL
  expect_message(ci <- confint(fe, parm = "delta"), "log(delta)",
                 fixed = TRUE)
  expect_identical(rownames(ci), "log(delta)")
  fa <- frm(y ~ x, family = acat(threshold = "equidistant"), data = d)
  expect_identical(lab(fa), c("x", "Intercept[1]", "delta"))
  fg <- frm(yg | thres(gr = g) ~ x, family = cumulative(), data = d)
  expect_identical(lab(fg)[c(2L, 6L)], c("Intercept[a,1]", "Intercept[b,1]"))
  # the old internal spelling still addresses the same row
  expect_identical(unname(confint(fg, parm = "tau_raw_5")),
                   unname(confint(fg, parm = "Intercept[b,1]")))
})

test_that("a threshold with no row of its own says where to read it", {
  # cumulative() holds its later thresholds as log increments, so
  # Intercept[2] is no internal parameter (review m7)
  d <- otn_data()
  fc <- frm(y ~ x, family = cumulative(), data = d)
  m <- tryCatch(confint(fc, parm = "Intercept[2]"),
                error = function(e) conditionMessage(e))
  expect_match(m, "fixef() reports every threshold", fixed = TRUE)
  # and the advice it gives answers
  h <- hypothesis(fc, "Intercept[2] = 0")
  expect_equal(h$hypothesis$Estimate,
               unname(fixef(fc)["Intercept[2]", "Estimate"]))
})

test_that("vcov_cluster() runs on an ordinal fit with disc held at 1", {
  # its template came from obj$env$parameters, where the map holds a
  # partly mapped component in a reduced shape: "A map factor length
  # must equal parameter length" on every ordinal fit without a modeled
  # disc, on 0.67.0 too (review m5, dev/fixes-rev-vcl.R)
  d <- otn_data()
  d$cl <- factor(rep(seq_len(20), length.out = nrow(d)))
  for (fam in list(cumulative(), sratio())) {
    fit <- frm(y ~ x, family = fam, data = d)
    V <- vcov_cluster(fit, cluster = d$cl, full = TRUE)
    expect_identical(rownames(V), rownames(vcov(fit, full = TRUE)))
    expect_true(all(is.finite(diag(V))) && all(diag(V) > 0))
    # the cluster scores add up to the gradient, zero at the optimum,
    # against the size of one cluster's score
    S <- cluster_scores(fit, d$cl)
    expect_lt(max(abs(colSums(S))), 1e-3 * max(abs(S)))
  }
})
