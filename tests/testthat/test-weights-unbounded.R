# A weighted likelihood with no maximum. When every row has a random
# effect of its own in mu, each density can narrow onto its own
# observation, and the likelihood grows like the peak height to the
# power sum(w) - n. The gated fuzz plan found it as Gamma + weights +
# ar1(tt + 0 | g) (one row per cell) + mo() + shape ~ z under
# profile = TRUE (seed 20502220), which died at "NA/NaN gradient
# evaluation"; without profile the same model reported a logLik of
# 1.4e11. mo(), the shape formula and profile play no part.

wu_data <- function(seed = 20261020) {
  set.seed(seed)
  ng <- 12; nt <- 6; n <- ng * nt
  d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                  tt = factor(rep(seq_len(nt), ng)),
                  x = rnorm(n), xs = rnorm(n), w = runif(n, 0.5, 2))
  d$obs <- factor(seq_len(n))
  eta <- 1 + 0.4 * d$x + rnorm(n, 0, 0.5)
  d$y <- rgamma(n, shape = 3, rate = 3 / exp(eta))
  d$ly <- log(d$y)
  d$k <- rpois(n, exp(eta))
  d$lsig <- rep(log(0.7), n)
  d$fac <- factor(rep(c("a", "b"), length.out = n))
  d
}

# the fit, with every warning it raised
wu_fit <- function(...) {
  warns <- character()
  fit <- withCallingHandlers(frm(...), warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(fit = fit, warns = warns)
}

wu_refusal <- "has no maximum"

test_that("an observation-level term with sum(w) > n is refused", {
  d <- wu_data()
  expect_gt(sum(d$w), nrow(d))
  G <- Gamma(link = "log")
  err <- tryCatch(frm(bf(y | weights(w) ~ x + ar1(tt + 0 | g)), family = G,
                      data = d), error = conditionMessage)
  expect_match(err, wu_refusal)
  # the term as written, and the condition for a maximum rather than a
  # factor to rescale the weights by
  expect_match(err, "ar1(tt + 0 | g)", fixed = TRUE)
  expect_match(err, "counts a row more than once", fixed = TRUE)
  expect_no_match(err, "Rescale", fixed = TRUE)
  # a prior that bounds the objective gives a degenerate estimate, so the
  # message names none as a remedy
  expect_no_match(err, "lognormal", fixed = TRUE)
  expect_match(err, "not the estimate informative", fixed = TRUE)
  expect_error(frm(bf(y | weights(w) ~ x + ar1(tt + 0 | g), shape ~ x),
                   family = G, data = d,
                   control = frmtmb_control(profile = TRUE)),
               wu_refusal)
  expect_error(frm(bf(y | weights(w) ~ x + (1 | obs)), family = G,
                   data = d, REML = TRUE),
               wu_refusal)
  # a gaussian response is the same model with sigma shrinking to 0;
  # the observation-level warning comes first
  allow_warnings(
    expect_error(frm(bf(ly | weights(w) ~ x + (1 | obs)),
                     family = gaussian(), data = d),
                 "sigma shrinks to 0"),
    "gives every observation its own random effect")
  # exgaussian narrows only with both of its scales
  allow_warnings(
    expect_error(frm(bf(ly | weights(w) ~ x + (1 | obs)),
                     family = exgaussian(), data = d),
                 "sigma and beta shrink to 0"),
    "gives every observation its own random effect")
})

test_that("only a prior that bounds the dispersion exempts the model", {
  d <- wu_data()
  fo <- bf(ly | weights(w) ~ x + (1 | obs))
  olre <- "gives every observation its own random effect"
  # finite at sigma = 0, or a slower power than the excess, or not on
  # sigma at all: the MAP objective still has no maximum
  for (p in list(set_prior("student_t(3, 0, 2.5)", class = "sigma"),
                 set_prior("gamma(2, 1)", class = "sigma"),
                 set_prior("normal(0, 5)", class = "b"))) {
    allow_warnings(
      expect_error(frm(fo, family = gaussian(), data = d, prior = p),
                   wu_refusal),
      olre)
  }
  # a density that vanishes at 0 faster than any power bounds it
  r <- allow_warnings(
    wu_fit(fo, family = gaussian(), data = d,
           prior = set_prior("lognormal(0, 1)", class = "sigma")),
    olre)
  expect_true(is.finite(as.numeric(logLik(r$fit))))
  # and so does a gamma tail on a shape
  r <- wu_fit(bf(y | weights(w) ~ x + (1 | obs)),
              family = Gamma(link = "log"), data = d,
              prior = set_prior("gamma(0.01, 0.01)", class = "shape"))
  expect_true(is.finite(as.numeric(logLik(r$fit))))
})

test_that("the refusal does not fire where a maximum exists", {
  d <- wu_data()
  G <- Gamma(link = "log")
  ok <- function(r) {
    expect_s3_class(r$fit, "frmtmb_fit")
    expect_true(is.finite(as.numeric(logLik(r$fit))))
    # a runaway fit stops on "false convergence" or a large gradient
    expect_identical(r$fit$opt$convergence, 0L)
    expect_length(r$warns, 0L)
  }
  # the weights sum to less than the number of rows
  dw <- d
  dw$w <- d$w * 0.7
  expect_lt(sum(dw$w), nrow(dw))
  rw <- wu_fit(bf(y | weights(w) ~ x + ar1(tt + 0 | g)), family = G,
               data = dw)
  ok(rw)
  # unit weights
  ru <- wu_fit(bf(y ~ x + ar1(tt + 0 | g)), family = G, data = d)
  ok(ru)
  # the weighted log shape stays within the unit-weight one's own size
  # of it (0.67 against 0.90 here); the runaway fits reached 33 to 172
  expect_lt(abs(rw$fit$estimates$betad[1] - ru$fit$estimates$betad[1]),
            abs(ru$fit$estimates$betad[1]))
  # no row owns its effect: six rows share each level
  ok(wu_fit(bf(y | weights(w) ~ x + (1 | g)), family = G, data = d))
  # a probability mass never exceeds one, whatever the weight
  fp <- frm(bf(k | weights(w) ~ x + (1 | obs)), family = poisson(),
            data = d)
  expect_true(is.finite(as.numeric(logLik(fp))))
  # the frame is not a fit, and the brms oracle reads it
  fr <- frm(bf(y | weights(w) ~ x + ar1(tt + 0 | g)), family = G,
            data = d, dry_run = "frame")
  expect_s3_class(fr, "frmtmb_frame")
})

test_that("a dispersion without an intercept cannot reach its limit", {
  d <- wu_data()
  olre <- "gives every observation its own random effect"
  # sigma ~ 0 + xs cannot shrink at every row at once, and an offset
  # alone has no free parameter: both likelihoods are bounded
  r1 <- allow_warnings(
    wu_fit(bf(ly | weights(w) ~ x + (1 | obs), sigma ~ 0 + xs),
           family = gaussian(), data = d),
    olre)
  expect_true(is.finite(as.numeric(logLik(r1$fit))))
  expect_identical(r1$fit$opt$convergence, 0L)
  r2 <- allow_warnings(
    wu_fit(bf(ly | weights(w) ~ x + (1 | obs), sigma ~ 0 + offset(lsig)),
           family = gaussian(), data = d),
    olre)
  r3 <- allow_warnings(
    wu_fit(bf(ly | weights(w) ~ x + (1 | obs), sigma = 0.7),
           family = gaussian(), data = d),
    olre)
  # the offset fit is the constant-sigma fit
  expect_equal(as.numeric(logLik(r2$fit)), as.numeric(logLik(r3$fit)))
  # cell means span the constant, so sigma ~ 0 + fac can reach the limit:
  # refused without a prior, and bounded by a proper prior on every
  # coefficient, as an intercept prior bounds it
  fo <- bf(ly | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac)
  allow_warnings(
    expect_error(frm(fo, family = gaussian(), data = d), wu_refusal),
    olre)
  r4 <- allow_warnings(
    wu_fit(fo, family = gaussian(), data = d,
           prior = set_prior("normal(0, 1)", class = "b", dpar = "sigma")),
    olre)
  expect_true(is.finite(as.numeric(logLik(r4$fit))))
  expect_identical(r4$fit$opt$convergence, 0L)
})

test_that("asym_laplace under Laplace has no maximum at any weight", {
  d <- wu_data()
  fo <- bf(ly ~ x + (1 | obs)) + asym_laplace()
  fr <- frm(fo, data = d, dry_run = "frame")
  # the Laplace approximation gives a row's own effect no curvature from
  # a log density that is linear on each side of its mode
  expect_error(check_weights_unbounded(fr$spec, fr, NULL, exact = FALSE),
               "Laplace approximation")
  # the exact integral is bounded at unit weights
  expect_null(check_weights_unbounded(fr$spec, fr, NULL, exact = TRUE))
  expect_error(frm(fo, data = d), "Laplace approximation")
})

test_that("the sampler's check asks whether the posterior integrates", {
  d <- wu_data()
  sig <- set_prior("student_t(3, 0, 2.5)", class = "sigma")
  frame_of <- function(dd) {
    allow_warnings(
      frm(bf(ly | weights(w) ~ x + (1 | obs)) + gaussian(), data = dd,
          dry_run = "frame"),
      "gives every observation its own random effect")
  }
  fr <- frame_of(d)
  # brms's default on sigma is finite at 0: improper once sum(w) - n >= 1
  expect_error(check_weights_unbounded(fr$spec, fr, sig, exact = TRUE,
                                       route = "sample"),
               "posterior of this model is improper")
  # with an excess below 1 the posterior integrates though the MAP
  # objective still has no maximum
  dh <- d
  dh$w <- 1
  dh$w[1] <- 1.5
  fh <- frame_of(dh)
  expect_null(check_weights_unbounded(fh$spec, fh, sig, exact = TRUE,
                                      route = "sample"))
  expect_error(check_weights_unbounded(fh$spec, fh, sig, exact = TRUE),
               wu_refusal)
  # a gamma tail on a shape makes it proper
  fg <- frm(bf(y | weights(w) ~ x + (1 | obs)) + Gamma(link = "log"),
            data = d, dry_run = "frame")
  shp <- set_prior("gamma(0.01, 0.01)", class = "shape")
  expect_null(check_weights_unbounded(fg$spec, fg, shp, exact = TRUE,
                                      route = "sample"))
})

test_that("cens() does not turn the check off: it counts uncensored rows", {
  # The check returned early on cens(), and the 0.73.0 families gained a
  # CDF, so a censoring column with no row censored switched the refusal
  # off (asym_laplace (1 | obs): logLik 762.92 at sigma 1.6e-05; Gamma
  # with weights 1.5: 5.1e7). A censored row's probability is at most 1
  # and its own effect can bring it to 1, so the uncensored rows still
  # run away.
  d <- wu_data()
  d$none <- "none"
  q <- stats::quantile(d$ly, 0.9)
  d$cc <- ifelse(d$ly > q, "right", "none")
  d$lyc <- pmin(d$ly, q)
  d$w15 <- 1.5
  expect_error(frm(bf(ly | cens(none) ~ x + (1 | obs)) + asym_laplace(),
                   data = d), "Laplace approximation")
  expect_gt(sum(d$cc == "right"), 0)
  expect_error(frm(bf(lyc | cens(cc) ~ x + (1 | obs)) + asym_laplace(),
                   data = d), "Laplace approximation")
  expect_error(frm(bf(y | weights(w15) + cens(none) ~ x + (1 | obs)) +
                     Gamma(link = "log"), data = d), wu_refusal)
  allow_warnings(
    expect_error(frm(bf(lyc | weights(w15) + cens(cc) ~ x + (1 | obs)) +
                       gaussian(), data = d),
                 "censored rows do not stop it"),
    "gives every observation its own random effect")
  # the absent case: six rows per level of g, not one, own no effect
  r <- wu_fit(bf(lyc | weights(w15) + cens(cc) ~ x + (1 | g)) + gaussian(),
              data = d)
  expect_true(is.finite(as.numeric(logLik(r$fit))))
  expect_identical(r$fit$opt$convergence, 0L)
})

test_that("zero_inflated_asym_laplace is asym_laplace on its nonzero rows", {
  # dev/rel073-zial-probe.R: with (1 | obs) the Laplace objective ran to
  # logLik 278.05 and 201.90 on 80 unweighted rows; a group of four rows
  # per level stayed at -119.29 and -124.25
  d <- wu_data()
  d$yz <- ifelse(seq_len(nrow(d)) %% 5 == 0, 0, d$ly)
  f <- zero_inflated_asym_laplace()
  expect_error(frm(bf(yz ~ x + (1 | obs)), family = f, data = d),
               "Laplace approximation")
  fr <- frm(bf(yz ~ x + (1 | obs)), family = f, data = d,
            dry_run = "frame")
  # the exact integral is bounded at unit weights
  expect_null(check_weights_unbounded(fr$spec, fr, NULL, exact = TRUE))
  r <- wu_fit(bf(yz ~ x + (1 | g)), family = f, data = d)
  expect_true(is.finite(as.numeric(logLik(r$fit))))
})
