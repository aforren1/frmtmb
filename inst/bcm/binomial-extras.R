# Binomial pieces the Bayesian Cognitive Modeling port needs and the
# core does not ship, written as ordinary rowwise custom families.
# vignette("bayesian-cognitive-modeling") shows them in place;
# tests/testthat/test-bcm-psychophysics.R,
# tests/testthat/test-bcm-signal-detection.R and
# tests/testthat/test-bcm-data-analysis.R fit this code. Load it with:
#
#   source(system.file("bcm", "binomial-extras.R", package = "frmtmb"))
#
# No family here needs frmtmb_structure(): every likelihood factorizes
# over rows, which is what a custom family is for. The structured
# protocol enters only where a row's contribution reads other rows; see
# inst/bcm/marginal.R.

## ---- bcm-cap1 ----
# min(1, u), on the tape.
#
# Several of the book's models write a probability as an expression that
# can exceed one and then cap it: the retention curve
# min(1, exp(-alpha t) + beta) and SIMPLE's min(1, sum(resp)). The cap
# is not decoration. Without it a cell where every item was recalled
# has a likelihood of theta^n with nothing holding theta down, so the
# maximum likelihood problem is UNBOUNDED and the fit walks off.
#
# RTMB has no pmin over advectors, and a branch on a parameter cannot be
# taped, but the identity
#
#   min(1, u) = u - max(0, u - 1) = u - (|u - 1| + (u - 1)) / 2
#
# is exact and taped, because abs() is. It has a kink at u = 1, where
# the derivative jumps from 1 to 0; that is a single point, it is where
# the model itself is kinked, and an optimizer approaching it from
# either side gets the derivative of the side it is on.
bcm_cap1 <- function(u) u - 0.5 * (abs(u - 1) + (u - 1))

## ---- bcm-probit-retired ----
# bcm_probit(), bcm_binomial_probit() and bcm_gaussian_probit() used to
# live here. The core link registry had no probit, so this file built
# one by hand and wrapped two families around it.
#
# The registry now holds probit, so all three are gone and the models
# they served are ordinary core families:
#
#   bcm_binomial_probit()  ->  binomial(link = "probit")
#   bcm_gaussian_probit()  ->  gaussian(link = "probit") with the known
#                              measurement SD on se(sd) instead of
#                              vreal(sd)
#
# Neither is a loss of faithfulness. The core binomial reads the
# registry's `logit_eta` for probit, which this file's hand-built link
# could not supply, so a saturated cell is now scored through
# dbinom_robust() rather than through a plain round trip. And a core
# gaussian is allowed se(), which is the channel the known standard
# deviation belonged in all along; the vreal() spelling was a way
# around the se() family gate that a custom family could not pass.

## ---- bcm-contaminant ----
# The contaminant binomial of Lee and Wagenmakers chapter 12.2. With
# probability phi a response comes from a contaminant process instead of
# the psychometric function, and the book gives that process a rate of
# its own with a beta(1, 1) prior, one per cell.
#
# That per-cell rate is integrated out here rather than estimated. Its
# only appearances are its own beta(1, 1) prior and one binomial
# likelihood, and a binomial rate with a beta(1, 1) prior integrates to
# the DISCRETE UNIFORM on 0..n:
#
#   integral over p of choose(n, r) p^r (1-p)^(n-r) dp = 1 / (n + 1)
#
# so the contaminant component contributes a flat -log(n + 1). Doing
# that integral is not a convenience. Left as a free parameter per cell,
# the contaminant fits every cell exactly, the mixture likelihood is
# maximized by phi = 1, and there is no maximum likelihood estimate of
# the psychometric function at all. The book's sampler never sees this
# because it integrates the same parameter out by averaging over it.
# The contamination rate takes a LOGIT link, where the book puts it on
# a probit. The reason is numerical and it is the seam the core's own
# binomial uses: `dpars[[".eta_phi"]]` is the linear predictor behind
# phi, so on a logit both log(phi) and log(1 - phi) can be recovered
# from it exactly, however far the predictor runs. Read off the natural
# scale instead, a probit rounds phi to 1 at a linear predictor of 8.3
# and log(1 - phi) becomes -Infinity where the true value is -8.3, which
# is what an optimizer walks into. The group distribution the book puts
# on the probit of phi becomes a group distribution on its logit; that
# is a different prior for the same model, and the test's Stan program
# says so.
bcm_contaminant <- function(link = "logit", link_phi = "logit") {
  frmtmb_family(
    "bcm_contaminant",
    dpars = c("mu", "phi"),
    links = list(mu = link, phi = link_phi),
    primary_dpars = "mu",
    type = "discrete",
    lpdf = function(y, dpars, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      phi <- dpars[["phi"]]
      eta <- dpars[[".eta_phi"]]
      if (is.null(eta)) {
        lphi <- log(phi)
        l1mphi <- log(1 - phi)
      } else {
        # -log(1 + exp(-eta)) and -log(1 + exp(eta)), exact at both ends
        lphi <- -RTMB::logspace_add(0 * eta, -eta)
        l1mphi <- -RTMB::logspace_add(0 * eta, eta)
      }
      # logspace_add keeps a decisive cell, where one component is many
      # nats above the other, from underflowing to a constant
      RTMB::logspace_add(
        l1mphi + RTMB::dbinom(y, size, dpars[["mu"]], log = TRUE),
        lphi - log(size + 1))
    },
    valid_y = function(y, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        stop("bcm_contaminant(): the response must be integer counts in ",
             "[0, trials]", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        size <- aterms[["trials"]]
        if (is.null(size)) size <- 1
        min(max(mean(y / size), 0.02), 0.98)
      },
      # a contamination rate that starts at a half cannot tell the two
      # components apart; the literature's own starting guess is small
      phi = function(y, aterms) 0.05
    ),
    post = list(
      mean_fn = function(dpars, aterms) {
        size <- aterms[["trials"]]
        if (is.null(size)) size <- 1
        phi <- dpars[["phi"]]
        # the contaminant's mean is the middle of 0..n
        size * ((1 - phi) * dpars[["mu"]] + phi / 2)
      },
      var_fn = function(dpars, aterms) {
        size <- aterms[["trials"]]
        if (is.null(size)) size <- 1
        phi <- dpars[["phi"]]
        mu <- dpars[["mu"]]
        m1 <- size * mu
        m2 <- size * mu * (1 - mu) + m1^2
        # the discrete uniform on 0..n has mean n/2 and variance
        # n (n + 2) / 12
        c1 <- size / 2
        c2 <- size * (size + 2) / 12 + c1^2
        mix <- (1 - phi) * m1 + phi * c1
        (1 - phi) * m2 + phi * c2 - mix^2
      }
    ),
    sim = function(dpars, aterms, n) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      size <- rep(as.numeric(size), length.out = n)
      cont <- stats::rbinom(n, 1L, dpars[["phi"]])
      ifelse(cont == 1,
             stats::rbinom(n, size, stats::runif(n)),
             stats::rbinom(n, size, dpars[["mu"]]))
    })
}

## ---- bcm-binomial-cdf ----
# The binomial with a CDF, which is all `cens()` asks a family for.
#
# ChaSaSoon (chapter 5.5) is the model that wants it: one observed score
# of 30 out of 50, and 949 earlier attempts known only to have scored
# somewhere between 15 and 25. That is interval censoring, and
# `cens("interval")` is the grammar's word for it.
#
# The core binomial has no `lcdf`, which is the one thing standing
# between it and a censored row, and the refusal says so. Four lines
# supply one, and the band is then written where a reader looks for it:
#
#   z | trials(n) + cens(cc, hi) + weights(w) ~ 1
#
# frmtmb reads a discrete censoring bound as INCLUSIVE, so the interval
# row scores F(25) - F(14) = P(15 <= Y <= 25), which is the book's Stan
# program exactly. `weights()` carries the 949 repeats.
#
# This family used to carry the band in its own density, because
# `cens()` refused every discrete family whatever CDF it supplied. That
# gate is gone; `dev/custom-findings.md` records what replaced it.
bcm_binomial_cdf <- function(link = "logit") {
  frmtmb_family(
    "bcm_binomial_cdf",
    dpars = "mu",
    links = list(mu = link),
    type = "discrete",
    accepts_aterms = c("trials", "weights", "cens", "trunc"),
    lpdf = function(y, dpars, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      RTMB::dbinom(y, size, dpars[["mu"]], log = TRUE)
    },
    lcdf = function(q, dpars, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      RTMB::pbinom(q, size, dpars[["mu"]])
    },
    valid_y = function(y, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        stop("bcm_binomial_cdf(): the response is an integer count in ",
             "[0, trials]", call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) {
        size <- aterms[["trials"]]
        if (is.null(size)) size <- 1
        min(max(mean(y / size), 0.02), 0.98)
      }
    ),
    post = list(
      mean_fn = function(dpars, aterms) {
        size <- aterms[["trials"]]
        if (is.null(size)) size <- 1
        size * dpars[["mu"]]
      },
      var_fn = function(dpars, aterms) {
        size <- aterms[["trials"]]
        if (is.null(size)) size <- 1
        size * dpars[["mu"]] * (1 - dpars[["mu"]])
      }
    ),
    sim = function(dpars, aterms, n) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      stats::rbinom(n, size, dpars[["mu"]])
    })
}

