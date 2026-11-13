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
  # with an excess below 1 the posterior integrates, and so does the MAP
  # objective have its maximum: frmtmb evaluates a prior on sigma on its
  # natural scale with the log link's Jacobian, sigma^1, so the
  # objective goes like sigma^(1 - excess) at 0. Through 0.81.0 the MAP
  # check left that power out and refused the model
  # (dev/reviews/wcheck-review-misc.R, part jac: slope -0.5 at E 0.5)
  dh <- d
  dh$w <- 1
  dh$w[1] <- 1.5
  fh <- frame_of(dh)
  expect_null(check_weights_unbounded(fh$spec, fh, sig, exact = TRUE,
                                      route = "sample"))
  expect_null(check_weights_unbounded(fh$spec, fh, sig, exact = TRUE))
  # at an excess above 1 both refuse
  dh$w[1:3] <- 1.6
  fh <- frame_of(dh)
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


test_that("a registered structure is checked by the rank it states", {
  # through 0.74.0 the check named the full-rank structures in a list,
  # so a structure another package registered was never one, and the
  # check failed open. It now reads the structure's full_rank. One that
  # states nothing leaves the check open: the covariance does not tell
  # whether a row's effect can reach any value (a latent with bounded
  # support keeps its maximum, below)
  local_frmtmb_registries()
  iid_lpdf <- function(b, theta, blk) {
    sum(RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE))
  }
  iid_vcov <- function(theta, blk) {
    matrix(exp(2 * theta[1]), 1, 1,
           dimnames = list(blk[["cnms"]], blk[["cnms"]]))
  }
  # an iid block of one coefficient under another name, full rank
  frmtmb_register_covstruct("myiid", lpdf = iid_lpdf, vcov = iid_vcov,
                            npar = 1L, sd_idx = 1L, full_rank = TRUE)
  # the same block, stating nothing
  frmtmb_register_covstruct("myiid2", lpdf = iid_lpdf, vcov = iid_vcov,
                            npar = 1L, sd_idx = 1L)
  # a block of two coefficients that moves only along (1, -1): its
  # covariance has rank one, and on columns a = b the effects add
  # nothing a row could use to reach its observation
  v <- c(1, -1) / sqrt(2)
  w <- c(1, 1) / sqrt(2)
  frmtmb_register_covstruct(
    "rankone",
    lpdf = function(b, theta, blk) {
      B <- RTMB::matrix(b, 2L)
      sum(RTMB::dnorm(as.vector(t(B) %*% v), 0, exp(theta[1]), log = TRUE)) +
        sum(RTMB::dnorm(as.vector(t(B) %*% w), 0, 1e-4, log = TRUE))
    },
    vcov = function(theta, blk) {
      V <- exp(2 * theta[1]) * tcrossprod(v)
      dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
      V
    },
    npar = 1L, sd_idx = 1L)
  # a latent with bounded support, a scaled beta(2, 2) on (-0.3, 0.3):
  # positive definite at every parameter value, and no runaway, since a
  # row's effect cannot reach its observation
  frmtmb_register_covstruct(
    "boxre",
    lpdf = function(b, theta, blk) {
      z <- (b + 0.3) / 0.6
      sum(log(6 * z * (1 - z) / 0.6))
    },
    vcov = function(theta, blk) {
      matrix(0.3^2 / 5, 1, 1, dimnames = list(blk[["cnms"]], blk[["cnms"]]))
    },
    npar = 0L)
  d <- wu_data()
  d$w15 <- 1.5
  d$one <- 1
  d$one2 <- 1
  G <- Gamma(link = "log")
  # the stated full-rank block is refused as the built-in one is
  e_us <- tryCatch(frm(bf(y | weights(w15) ~ x + (1 | obs)), family = G,
                       data = d), error = conditionMessage)
  e_my <- tryCatch(frm(bf(y | weights(w15) ~ x + myiid(1 | obs)),
                       family = G, data = d), error = conditionMessage)
  expect_match(e_my, wu_refusal)
  expect_identical(sub("myiid(1 | obs)", "(1 | obs)", e_my, fixed = TRUE),
                   e_us)
  expect_true(frmtmb_covstruct_info("myiid")$full_rank)
  # the unstated one is not refused: the check is open for it, and the
  # reader says so
  expect_identical(frmtmb_covstruct_info("myiid2")$full_rank, NA)
  expect_false(frmtmb:::cs_full_rank(list(covstruct = "myiid2", dim = 1L,
                                          cnms = "Intercept",
                                          theta_idx = 1L)))
  # the bounded latent, unstated, fits and has its maximum
  rb <- wu_fit(bf(y | weights(w15) ~ x + boxre(1 | obs)), family = G,
               data = d)
  expect_true(is.finite(as.numeric(logLik(rb$fit))))
  expect_lt(exp(rb$fit$estimates$betad[1L]), 100)
  # the rank-one block, unstated, is not refused, and the model has its
  # maximum
  r <- wu_fit(bf(y | weights(w15) ~ x + rankone(0 + one + one2 | obs)),
              family = G, data = d)
  # the block adds nothing a row can use, so the maximum is the one
  # without it, to the 1e-4 sd its unused direction keeps; its own sd
  # is flat there, which the fit says in its warnings (wu_fit() keeps
  # them)
  r0 <- wu_fit(bf(y | weights(w15) ~ x), family = G, data = d)
  ll <- as.numeric(logLik(r$fit))
  expect_true(is.finite(ll))
  expect_equal(ll, as.numeric(logLik(r0$fit)), tolerance = 1e-6)
  expect_equal(unname(fixef(r$fit)[, 1]), unname(fixef(r0$fit)[, 1]),
               tolerance = 1e-4)
  expect_true(frmtmb_covstruct_info("us")$full_rank)
  expect_false(frmtmb_covstruct_info("rr")$full_rank)
  # a name nothing registered: the reader says so, in its own words
  expect_error(frmtmb_covstruct_info("nosuch"),
               "no covariance structure 'nosuch' is registered in this session",
               fixed = TRUE, class = "frmtmb_error")
  # a slot that is a function of `...` takes any arguments
  frmtmb_register_covstruct("dotsiid", lpdf = function(...) iid_lpdf(...),
                            vcov = function(...) iid_vcov(...),
                            npar = 1L, sd_idx = 1L)
  expect_identical(frmtmb_covstruct_info("dotsiid")$full_rank, NA)
})

test_that("a Student-t latent per row is full rank, and refused under
           weights that sum to more than the rows", {
  # through 0.74.0 us_t and diag_t were not full rank, so the check was
  # open for them: a weighted Gamma model with (1 | gr(obs, dist =
  # "student")) ran to a logLik of 4.7e7, the shape to 1e7. A t density
  # is positive everywhere, so a row's effect reaches its observation as
  # a gaussian one does
  d <- wu_data()
  d$w15 <- 1.5
  G <- Gamma(link = "log")
  expect_true(frmtmb_covstruct_info("us_t")$full_rank)
  expect_true(frmtmb_covstruct_info("diag_t")$full_rank)
  err <- tryCatch(frm(y | weights(w15) ~ x + (1 | gr(obs, dist = "student")),
                      family = G, data = d), error = conditionMessage)
  expect_match(err, wu_refusal)
  expect_match(err, "gr(obs, dist = \"student\")", fixed = TRUE)
  # the absent case: unweighted, the same model has its maximum
  r <- wu_fit(y ~ x + (1 | gr(obs, dist = "student")), family = G, data = d)
  expect_true(is.finite(as.numeric(logLik(r$fit))))
  expect_lt(exp(r$fit$estimates$betad[1L]), 100)
})

# the fit, with every warning and message it raised
wu_fit_all <- function(...) {
  msgs <- character()
  r <- withCallingHandlers(wu_fit(...), message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  r$msgs <- msgs
  r
}

test_that("each response of a multivariate model is checked on its rows", {
  # through 0.81.0 the check returned for more than one response, and
  # mvbf() with (1 | obs) per response ran away without a refusal: two
  # gaussian responses reached logLik 836.78 at log sigma -37.49 and
  # -21.39 with false convergence (dev/wcheck-cases.R, part mv)
  d <- wu_data()
  set.seed(20261111)
  d$ly2 <- d$ly + stats::rnorm(nrow(d), 0, 0.5)
  d$k2 <- stats::rpois(nrow(d), exp(0.3 + 0.2 * d$x +
                                     stats::rnorm(nrow(d), 0, 0.6)))
  olre <- "gives every observation its own random effect"
  G <- gaussian()
  err <- allow_warnings(
    tryCatch(frm(mvbf(bf(ly ~ x + (1 | obs)),
                      bf(ly2 | weights(w) ~ x + (1 | obs))),
                 family = G, data = d), error = conditionMessage),
    olre)
  expect_match(err, wu_refusal)
  # the response that runs away, and its term as its formula writes it
  expect_match(err, "in the mu of response ly2", fixed = TRUE)
  expect_match(err, "`(1 | obs)` gives", fixed = TRUE)
  # a count response beside it does not hide a Gamma one
  expect_error(frm(mvbf(bf(k ~ x + (1 | obs)) + poisson(),
                        bf(y | weights(w) ~ x + (1 | obs)) +
                          Gamma(link = "log")), data = d),
               "in the mu of response y")
  # nor does a block that the responses share through an ID
  allow_warnings(
    expect_error(frm(mvbf(bf(ly | weights(w) ~ x + (1 | p | obs)),
                          bf(ly2 | weights(w) ~ x + (1 | p | obs))),
                     family = G, data = d), wu_refusal),
    olre)
  # the exact integral of quadrature = TRUE has no maximum either
  allow_warnings(
    expect_error(frm(mvbf(bf(ly | weights(w) ~ x + (1 | obs)),
                          bf(ly2 | weights(w) ~ x + (1 | obs))),
                     family = G, data = d, quadrature = TRUE), wu_refusal),
    olre)
  # the absent case: an overdispersed count with an effect per row is
  # identified, weighted or not, and fits with nothing to say
  mv_counts <- list(
    mvbf(bf(k ~ x + (1 | obs)), bf(k2 ~ x + (1 | obs))),
    mvbf(bf(k | weights(w) ~ x + (1 | obs)),
         bf(k2 | weights(w) ~ x + (1 | obs))))
  for (fo in mv_counts) {
    r <- wu_fit_all(fo, family = poisson(), data = d)
    expect_true(is.finite(as.numeric(logLik(r$fit))))
    expect_identical(r$fit$opt$convergence, 0L)
    expect_length(r$warns, 0L)
    expect_length(r$msgs, 0L)
  }
  # and gaussian responses whose weights sum to fewer than the rows
  dw <- d
  dw$w <- d$w * 0.7
  expect_lt(sum(dw$w), nrow(dw))
  r <- wu_fit_all(mvbf(bf(ly | weights(w) ~ x + (1 | obs)),
                       bf(ly2 | weights(w) ~ x + (1 | obs))),
                  family = G, data = dw)
  expect_true(is.finite(as.numeric(logLik(r$fit))))
  expect_identical(r$fit$opt$convergence, 0L)
  expect_true(all(grepl(olre, r$warns, fixed = TRUE)))
  expect_length(r$msgs, 0L)
})

test_that("a censored row costs nothing under Laplace when a mean brings
           every one to 1", {
  # 0.81.0 took 1/2 off the power per censored row, so an uncensored
  # excess below half the censored rows was let through: on 80 rows with
  # 16 right censored, excesses of 0.5 to 4 stopped at an interior
  # estimate and 7.5 and 8 ran away (log sigma -14.16, -13.64, conv 1).
  # A mean above every bound puts each censored row's probability at 1
  # at its own effect's prior mean, where the data give that effect no
  # curvature, and the objective grows with the uncensored excess alone
  # (dev/wcheck-cens-path.R: slopes 0.5, 2, 4 and 7.5)
  d <- wu_data()
  olre <- "gives every observation its own random effect"
  q <- unname(stats::quantile(d$ly, 0.8))
  d$cc <- ifelse(d$ly > q, "right", "none")
  d$lyc <- pmin(d$ly, q)
  unc <- d$cc == "none"
  E <- 2
  d$we <- ifelse(unc, 1 + E / sum(unc), 1)
  # below the 0.81.0 allowance
  expect_lt(E, sum(!unc) / 2)
  fo <- bf(lyc | weights(we) + cens(cc) ~ x + (1 | obs))
  err <- allow_warnings(
    tryCatch(frm(fo, family = gaussian(), data = d),
             error = conditionMessage),
    olre)
  expect_match(err, wu_refusal)
  expect_match(err, "censored rows do not stop it", fixed = TRUE)
  expect_no_match(err, "take 1/2 each", fixed = TRUE)
  # the growth itself, along that mean, with the effects' sd held
  slope <- function(dd) {
    ob <- allow_warnings(frm(fo, family = gaussian(), data = dd,
                             dry_run = "objective"), olre)
    ob <- if (!is.null(ob$obj)) ob$obj else ob
    p <- ob$par
    p[names(p) == "beta"] <- c(q + 2, 0)
    p[names(p) == "theta"] <- log(stats::sd(dd$lyc - q - 2))
    at <- function(ls) {
      p[names(p) == "betad"] <- ls
      ob$fn(p)
    }
    (at(-10) - at(-20)) / 10
  }
  expect_lt(abs(slope(d) - E), 1e-3 * E)
  # the control: unit weights have no excess, and the same path is flat
  d0 <- d
  d0$we <- 1
  expect_lt(abs(slope(d0)), 1e-3 * E)
  r0 <- allow_warnings(wu_fit(fo, family = gaussian(), data = d0), olre)
  expect_true(is.finite(as.numeric(logLik(r0$fit))))
  # asym_laplace with all rows but 8 right censored: 8 rows grow, and
  # 0.81.0 took 32 off
  da <- d
  o <- order(da$ly, decreasing = TRUE)
  k <- nrow(da) - 8L
  da$ca <- "none"
  da$ca[o[seq_len(k)]] <- "right"
  da$lyac <- pmin(da$ly, da$ly[o[k]])
  expect_error(frm(lyac | cens(ca) ~ x + (1 | obs), family = asym_laplace(),
                   data = da), "Laplace approximation")
  # both sides censored: one mean cannot bring both to 1, so the smaller
  # side still costs 1/2 per row, and the message counts it
  ql <- unname(stats::quantile(d$ly, 0.1))
  d$c2 <- ifelse(d$ly > q, "right", ifelse(d$ly < ql, "left", "none"))
  d$ly2c <- pmin(pmax(d$ly, ql), q)
  n_left <- sum(d$c2 == "left")
  unc2 <- d$c2 == "none"
  d$w2 <- ifelse(unc2, 1 + (n_left / 2 + 1) / sum(unc2), 1)
  err2 <- allow_warnings(
    tryCatch(frm(bf(ly2c | weights(w2) + cens(c2) ~ x + (1 | obs)),
                 family = gaussian(), data = d), error = conditionMessage),
    olre)
  expect_match(err2, paste0("the ", n_left, " that a constant mean cannot"),
               fixed = TRUE)
})

test_that("a free coefficient that moves one cell's rows alone is refused", {
  # through 0.81.0 a prior that bounded some cell means of sigma ~ 0 +
  # fac, or the intercept of sigma ~ 1 + fac, exempted the model while
  # the free cell ran away (logLik 377.18 at log sigma -36.89; 437.86 at
  # -23.24; dev/wcheck-cases.R, part prior)
  d <- wu_data()
  olre <- "gives every observation its own random effect"
  expect_gt(sum(d$w[d$fac == "a"] - 1), 0)
  expect_gt(sum(d$w[d$fac == "b"] - 1), 0)
  P <- function(...) set_prior(..., dpar = "sigma")
  fo <- bf(ly | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac)
  err <- allow_warnings(
    tryCatch(frm(fo, family = gaussian(), data = d,
                 prior = P("normal(0, 1)", class = "b", coef = "faca")),
             error = conditionMessage),
    olre)
  expect_match(err, wu_refusal)
  expect_match(err, "the coefficient `facb` of sigma moves alone",
               fixed = TRUE)
  # a student_t tail on every cell bounds none of them
  allow_warnings(
    expect_error(frm(fo, family = gaussian(), data = d,
                     prior = P("student_t(3, 0, 2.5)", class = "b")),
                 wu_refusal),
    olre)
  # sigma ~ 1 + fac: a free treatment contrast alone runs away when its
  # cell's weights sum above its rows. With a prior on the Intercept it
  # need not: frmtmb puts that prior on the centered intercept, which
  # the contrast moves too, so the check stays open there. The model
  # below has its maximum with the prior (logLik -94.2889 on 0.81.0; a
  # profile in b_facb falls to -193.36 at -35,
  # dev/reviews/wcheck-review-center3.R) and none without it (slope 2 =
  # the cell's excess)
  set.seed(20261111)
  dc <- data.frame(obs = factor(seq_len(80)), x = stats::rnorm(80),
                   fac = factor(rep(c("a", "b"), 40)))
  dc$y <- 0.5 + 0.4 * dc$x +
    stats::rnorm(80, 0, ifelse(dc$fac == "a", 0.6, 1.2))
  dc$w <- ifelse(dc$fac == "a", 0.8, 1.05)
  fc <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
  allow_warnings(
    expect_error(frm(fc, family = gaussian(), data = dc),
                 "the coefficient `facb` of sigma moves alone"),
    olre)
  r1 <- wu_fit(fc, family = gaussian(), data = dc,
               prior = P("normal(0, 1)", class = "Intercept"))
  expect_true(is.finite(as.numeric(logLik(r1$fit))))
  expect_identical(r1$fit$opt$convergence, 0L)
  expect_true(all(grepl(olre, r1$warns, fixed = TRUE)))
  # the absent case: every cell bounded, and a slope that reaches no
  # row's limit alone
  r <- allow_warnings(
    wu_fit(bf(ly | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac + xs),
           family = gaussian(), data = d,
           prior = c(P("normal(0, 1)", class = "b", coef = "faca"),
                     P("normal(0, 1)", class = "b", coef = "facb"))),
    olre)
  expect_true(is.finite(as.numeric(logLik(r$fit))))
  expect_identical(r$fit$opt$convergence, 0L)
})

test_that("a power-law prior tail on a shape bounds only past the growth", {
  # a Gamma density narrows like shape^(-1/2), so the likelihood grows
  # like shape^(excess / 2); a student_t(nu) prior falls like
  # shape^-(nu + 1), and its Jacobian on the log scale adds shape^1, so
  # the objective goes like shape^(excess / 2 - nu): bounded when nu
  # exceeds excess / 2. Through 0.81.0 any student_t on the shape
  # exempted the model, and student_t(3, 0, 2.5) ran to logLik 1.8e10 at
  # log shape 51.7 (dev/wcheck-cases.R, part prior)
  d <- wu_data()
  E <- sum(d$w) - nrow(d)
  expect_gt(E / 2, 4)
  expect_lt(E / 2, 30)
  fo <- bf(y | weights(w) ~ x + (1 | obs))
  G <- Gamma(link = "log")
  sh <- function(p) set_prior(p, class = "shape")
  expect_error(frm(fo, family = G, data = d,
                   prior = sh("student_t(3, 0, 2.5)")), wu_refusal)
  expect_error(frm(fo, family = G, data = d, prior = sh("inv_gamma(2, 1)")),
               wu_refusal)
  # half a power on each side of the boundary nu = excess / 2: below it
  # the objective grows like shape^0.5 and is refused, above it it falls
  # like shape^-0.5 and is not (dev/wcheck-fr-boundary.R: slopes +0.500
  # and -0.500 at E 12 and nu 5.5, 6.5)
  t_at <- function(nu) sh(sprintf("student_t(%.6f, 0, 2.5)", nu))
  expect_error(frm(fo, family = G, data = d, prior = t_at(E / 2 - 0.5)),
               wu_refusal)
  ra <- wu_fit(fo, family = G, data = d, prior = t_at(E / 2 + 0.5))
  expect_s3_class(ra$fit, "frmtmb_fit")
  # the absent case: a tail faster than the growth
  r <- wu_fit(fo, family = G, data = d, prior = sh("student_t(30, 0, 2.5)"))
  expect_true(is.finite(as.numeric(logLik(r$fit))))
  expect_identical(r$fit$opt$convergence, 0L)
  # the posterior integrates when nu exceeds the growth's power
  fg <- frm(fo, family = G, data = d, dry_run = "frame")
  expect_null(check_weights_unbounded(fg$spec, fg,
                                      sh("student_t(30, 0, 2.5)"),
                                      exact = TRUE, route = "sample"))
  expect_error(check_weights_unbounded(fg$spec, fg,
                                       sh("student_t(3, 0, 2.5)"),
                                       exact = TRUE, route = "sample"),
               "posterior of this model is improper")
})

test_that("a prior on sigma counts its Jacobian, as the objective does", {
  # frmtmb evaluates a class sigma prior on the natural scale with the
  # log link's Jacobian, so the MAP objective goes like
  # sigma^(order + 1 - excess) at 0: bounded when order + 1 exceeds the
  # excess. Through 0.81.0 the check read order > excess and refused
  # bounded models (dev/reviews/wcheck-review-misc.R, part jac: slope
  # -0.5 per unit of -log sigma at an excess of 0.5 under normal(0, 1))
  d <- wu_data()
  n <- nrow(d)
  olre <- "gives every observation its own random effect"
  fo <- bf(ly | weights(we) ~ x + (1 | obs))
  sg <- function(p) set_prior(p, class = "sigma")
  at <- function(E, p) {
    d$we <- 1 + E / n
    allow_warnings(
      tryCatch(wu_fit(fo, family = gaussian(), data = d, prior = sg(p))$fit,
               error = conditionMessage),
      olre)
  }
  # order 0: bounded below an excess of 1, not above it
  for (p in c("normal(0, 1)", "student_t(3, 0, 2.5)", "exponential(1)")) {
    expect_s3_class(at(0.5, p), "frmtmb_fit")
    expect_match(at(1.5, p), wu_refusal)
  }
  # gamma(a, 1) has order a - 1: bounded when a exceeds the excess
  expect_s3_class(at(2, "gamma(2.5, 1)"), "frmtmb_fit")
  expect_match(at(2, "gamma(1.5, 1)"), wu_refusal)
})

test_that("a prior that bounds the mean keeps censored rows' cost", {
  # a constant mean puts every right-censored row's probability at 1
  # only if it can rise above the bound; uniform(-5, 0.5) on the
  # intercept keeps it below 1.668, so each censored row still takes
  # 1/2 and the model has its maximum (dev/reviews/wcheck-review-cens2.R:
  # profile slopes -6.33 and -2.33 at excesses 2 and 6)
  set.seed(20261111)
  d <- data.frame(obs = factor(seq_len(80)), x = stats::rnorm(80))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(80)
  q8 <- unname(stats::quantile(d$y, 0.8))
  d$cc <- ifelse(d$y > q8, "right", "none")
  d$yc <- pmin(d$y, q8)
  nu <- sum(d$cc == "none")
  olre <- "gives every observation its own random effect"
  fo <- bf(yc | weights(w) + cens(cc) ~ 1 + (1 | obs))
  pr <- set_prior("uniform(-5, 0.5)", class = "Intercept", lb = -5,
                  ub = 0.5)
  for (E in c(2, 6)) {
    d$w <- ifelse(d$cc == "none", 1 + E / nu, 1)
    r <- wu_fit(fo, family = gaussian(), data = d, prior = pr)
    expect_true(is.finite(as.numeric(logLik(r$fit))))
    expect_identical(r$fit$opt$convergence, 0L)
    # without the prior the mean rises above the bound, and it runs away
    allow_warnings(
      expect_error(frm(fo, family = gaussian(), data = d), wu_refusal),
      olre)
  }
})

test_that("a prior names a response in brms's spelling too", {
  # resp = "ya" is brms's name for the response y_a, and frmtmb applies
  # it; the check read only "y_a" and refused a model that a prior on
  # each sigma bounds (dev/reviews/wcheck-review-respname.R)
  set.seed(20261021)
  n <- 80
  d <- data.frame(obs = factor(seq_len(n)), x = stats::rnorm(n),
                  w = stats::runif(n, 0.5, 2))
  d$y_a <- 0.5 + 0.4 * d$x + stats::rnorm(n)
  d$y.b <- d$y_a + stats::rnorm(n)
  expect_gt(sum(d$w), n)
  fo <- mvbf(bf(y_a | weights(w) ~ x + (1 | obs)),
             bf(y.b | weights(w) ~ x + (1 | obs)))
  olre <- "gives every observation its own random effect"
  for (rs in list(c("ya", "yb"), c("y_a", "y.b"))) {
    pr <- c(set_prior("lognormal(0, 1)", class = "sigma", resp = rs[1]),
            set_prior("lognormal(0, 1)", class = "sigma", resp = rs[2]))
    r <- wu_fit(fo, family = gaussian(), data = d, prior = pr)
    expect_true(is.finite(as.numeric(logLik(r$fit))))
  }
  # and with no prior the first response is refused by its own name
  allow_warnings(
    expect_error(frm(fo, family = gaussian(), data = d),
                 "in the mu of response y_a"),
    olre)
})

test_that("the Jacobian a natural-scale prior carries follows the link", {
  # log |dx / d eta| goes like x^j: j = 1 for log at both ends; 0 for
  # identity; softplus 1 toward 0 and 0 toward infinity; squareplus 2
  # toward 0 and 0 toward infinity. A first fix applied the log link's
  # j = 1 to every link and refused bounded models
  # (dev/reviews/wcheck-short-shapelinks.R, wcheck-short-links.R;
  # every boundary profiled on both sides in dev/wcheck-fr-links.R)
  set.seed(20261111)
  n <- 80
  d <- data.frame(obs = factor(seq_len(n)), x = stats::rnorm(n))
  mu <- exp(0.5 + 0.3 * d$x + stats::rnorm(n, 0, 0.3))
  d$yg <- stats::rgamma(n, shape = 4, rate = 4 / mu)
  d$y <- 0.5 + 0.4 * d$x + stats::rnorm(n)
  olre <- "gives every observation its own random effect"
  # Gamma, excess 12 (growth shape^6): student_t(5.5) falls like
  # shape^-6.5, so the objective is bounded under a link with j = 0
  # (0.81.0: logLik -95.1383, conv 0) and not under log (j = 1)
  d$w <- 1 + 12 / n
  fg <- bf(yg | weights(w) ~ x + (1 | obs))
  pt <- set_prior("student_t(5.5, 0, 2.5)", class = "shape")
  for (ln in c("identity", "softplus")) {
    r <- wu_fit(fg, family = brmsfamily("Gamma", link = "log",
                                        link_shape = ln),
                data = d, prior = pt)
    expect_true(is.finite(as.numeric(logLik(r$fit))))
    expect_identical(r$fit$opt$convergence, 0L)
  }
  expect_error(frm(fg, family = brmsfamily("Gamma", link = "log",
                                           link_shape = "log"),
                   data = d, prior = pt), wu_refusal)
  # gaussian, normal(0, 1) on sigma (order 0): squareplus's j = 2 keeps
  # an excess of 1.5 bounded (slope -0.5 per unit of -log sigma), where
  # 0.81.0 refused it; identity's j = 0 leaves an excess of 0.5 unbounded
  # (slope +0.5)
  fn <- bf(y | weights(w) ~ x + (1 | obs))
  ps <- set_prior("normal(0, 1)", class = "sigma")
  ds <- d
  ds$w <- 1 + 1.5 / n
  r <- allow_warnings(
    wu_fit(fn, family = brmsfamily("gaussian", link_sigma = "squareplus"),
           data = ds, prior = ps),
    olre)
  expect_true(is.finite(as.numeric(logLik(r$fit))))
  expect_identical(r$fit$opt$convergence, 0L)
  di <- d
  di$w <- 1 + 0.5 / n
  allow_warnings(
    expect_error(frm(fn, family = brmsfamily("gaussian",
                                             link_sigma = "identity"),
                     data = di, prior = ps), wu_refusal),
    olre)
  # frechet's nu, default link logm1 (dnu/deta = nu - 1, j = 1 as nu
  # grows), excess 6: student_t(5.5) on nu leaves the MAP objective
  # growing like nu^0.5, student_t(6.5) bounds it
  # (dev/wcheck-fr3-frechet.R: slopes +0.50 and -0.50)
  set.seed(20261111)
  df_ <- data.frame(obs = factor(seq_len(n)), x = stats::rnorm(n))
  mf <- exp(0.5 + 0.3 * df_$x + stats::rnorm(n, 0, 0.3))
  df_$yf <- mf * (-log(stats::runif(n)))^(-1 / 4)
  df_$w <- 1 + 6 / n
  ff <- frm(bf(yf | weights(w) ~ x + (1 | obs)), family = frechet(),
            data = df_, dry_run = "frame")
  pn <- function(nu) {
    set_prior(sprintf("student_t(%g, 0, 2.5)", nu), class = "nu")
  }
  expect_error(check_weights_unbounded(ff$spec, ff, pn(5.5)), wu_refusal)
  expect_null(check_weights_unbounded(ff$spec, ff, pn(6.5)))
})

test_that("the sampler's verdict does not depend on the link", {
  # with the Jacobian the posterior's integral over eta is the integral
  # of p(x) L(x) dx, so the sample route reads power 1 for every link,
  # as 0.81.0 did; only the MAP objective carries the link's power
  # (dev/reviews/wcheck-short2-routes.R; the profiled posterior mass of
  # the first case falls 0.50 per decade of sigma,
  # dev/reviews/wcheck-short2-mass.R)
  set.seed(20261111)
  n <- 80
  d <- data.frame(obs = factor(seq_len(n)), x = stats::rnorm(n))
  d$y <- 0.5 + 0.4 * d$x + stats::rnorm(n)
  mu <- exp(0.5 + 0.3 * d$x + stats::rnorm(n, 0, 0.3))
  d$yg <- stats::rgamma(n, shape = 4, rate = 4 / mu)
  olre <- "gives every observation its own random effect"
  improper <- "posterior of this model is improper"
  frame_of <- function(fo, fam, w) {
    d$w <- w
    allow_warnings(frm(fo, family = fam, data = d, dry_run = "frame"),
                   olre)
  }
  fn <- bf(y | weights(w) ~ x + (1 | obs))
  ps <- set_prior("normal(0, 1)", class = "sigma")
  # identity on sigma, excess 0.5: proper (the MAP objective is not
  # bounded, and the MAP route refuses it)
  fi <- frame_of(fn, brmsfamily("gaussian", link_sigma = "identity"),
                 1 + 0.5 / n)
  expect_null(check_weights_unbounded(fi$spec, fi, ps, exact = TRUE,
                                      route = "sample"))
  expect_error(check_weights_unbounded(fi$spec, fi, ps, exact = TRUE),
               wu_refusal)
  # squareplus on sigma, excess 1.5: improper (the MAP objective is
  # bounded, and the MAP route passes it)
  fs <- frame_of(fn, brmsfamily("gaussian", link_sigma = "squareplus"),
                 1 + 1.5 / n)
  expect_error(check_weights_unbounded(fs$spec, fs, ps, exact = TRUE,
                                       route = "sample"), improper)
  expect_null(check_weights_unbounded(fs$spec, fs, ps, exact = TRUE))
  # Gamma, excess 12, student_t(5.5) on the shape: improper under every
  # link
  fg <- bf(yg | weights(w) ~ x + (1 | obs))
  pt <- set_prior("student_t(5.5, 0, 2.5)", class = "shape")
  for (ln in c("identity", "softplus", "squareplus")) {
    fr <- frame_of(fg, brmsfamily("Gamma", link = "log", link_shape = ln),
                   1 + 12 / n)
    expect_error(check_weights_unbounded(fr$spec, fr, pt, exact = TRUE,
                                         route = "sample"), improper)
  }
})
