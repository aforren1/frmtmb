# Binomial pieces the Bayesian Cognitive Modeling port needs and the
# core does not ship, written as an ordinary link object and ordinary
# rowwise custom families.
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

## ---- bcm-probit ----
# The probit link.
#
# frmtmb's link registry has identity, log, logit, cloglog, inverse,
# logm1, tan_half and power12, and no probit; brms has one. Signal
# detection theory is DEFINED on the normal scale, and so are the latent
# traits of the multinomial processing trees, so the port needs it.
#
# A link is a plain list of four functions, so supplying one costs
# nothing and needs no change to the core. The RTMB spellings are the
# taped ones: stats::pnorm clamps at C level in ways the tape cannot
# see, which is the reason the core keeps its own registry rather than
# calling stats::make.link.
bcm_probit <- function() {
  list(
    name = "probit",
    linkfun = function(mu) RTMB::qnorm(mu),
    linkinv = function(eta) RTMB::pnorm(eta),
    mu_eta = function(eta) RTMB::dnorm(eta))
}

## ---- bcm-binomial-probit ----
# A binomial on the probit scale, which is what an equal-variance
# Gaussian signal detection model is. Identical to the core binomial
# except for the link, and without the core's `logit_eta` shortcut: a
# probit has no exact log-odds form to recover from its linear
# predictor, so a saturated cell is scored through the plain round trip
# and a design that saturates has to be caught by the fit, not hidden.
bcm_binomial_probit <- function() {
  frmtmb_family(
    "bcm_binomial_probit",
    dpars = "mu",
    links = list(mu = bcm_probit()),
    type = "discrete",
    lpdf = function(y, dpars, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      RTMB::dbinom(y, size, dpars[["mu"]], log = TRUE)
    },
    valid_y = function(y, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        stop("bcm_binomial_probit(): the response must be integer counts ",
             "in [0, trials]", call. = FALSE)
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

## ---- bcm-binomial-band ----
# A binomial whose response is a BAND of counts rather than one count.
#
# ChaSaSoon (chapter 5.5) is the model that wants it: one observed score
# of 30 out of 50, and 949 earlier attempts known only to have scored
# somewhere between 15 and 25. That is interval censoring, and
# `cens("interval")` is the grammar's word for it, but it does not reach
# here. TWO gates stand in the way and only the first is a seam:
#
#   1. `cens()` needs a family with a CDF, and the core binomial has
#      none. That gate says how to pass it, and `lcdf` below does.
#   2. `cens()` refuses ANY family whose `type` is "discrete", CDF or
#      not: "cens() is not supported for discrete families yet
#      (truncation is)". A custom family cannot opt in.
#
# So the band is the family's own business instead. The response is the
# SMALLEST count consistent with the row and `vint(hi)` the largest, so
#
#   log P(y <= Y <= hi) = log(F(hi) - F(y - 1))
#
# which is the exact density when hi == y and the band probability
# otherwise, with no branch at all. `weights()` then carries how many
# attempts fell in the band. The gate is recorded in
# dev/bcm-findings.md.
#
# `lcdf` is kept even though nothing in the core reaches it: it is the
# same arithmetic the band uses, and a family that supplies it is what
# the first gate asks for.
bcm_binomial_band <- function(link = "logit") {
  frmtmb_family(
    "bcm_binomial_band",
    dpars = "mu",
    links = list(mu = link),
    type = "discrete",
    required_aterms = "vint1",
    lpdf = function(y, dpars, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      log(RTMB::pbinom(aterms[["vint1"]], size, dpars[["mu"]]) -
            RTMB::pbinom(y - 1, size, dpars[["mu"]]))
    },
    lcdf = function(q, dpars, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      RTMB::pbinom(q, size, dpars[["mu"]])
    },
    valid_y = function(y, aterms) {
      size <- aterms[["trials"]]
      if (is.null(size)) size <- 1
      hi <- aterms[["vint1"]]
      if (any(y < 0) || any(y > size) || any(y != round(y))) {
        stop("bcm_binomial_band(): the response is the smallest count ",
             "consistent with the row, an integer in [0, trials]",
             call. = FALSE)
      }
      if (any(hi < y) || any(hi > size) || any(hi != round(hi))) {
        stop("bcm_binomial_band(): vint(hi) is the largest count ",
             "consistent with the row, so it is an integer in ",
             "[response, trials]", call. = FALSE)
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

## ---- bcm-gaussian-probit ----
# A gaussian whose mean is Phi of its linear predictor, with a KNOWN
# standard deviation carried per row.
#
# The ESP replication of chapter 16.2 needs it. There a person's
# extraversion score is normal around 100 * Phi(theta) with a known
# measurement standard deviation, and theta is one half of a correlated
# latent pair whose other half is a probit rate. Dividing the response
# by 100 makes the mean Phi(theta) exactly, so the whole nonlinearity is
# a link and the latent pair is an ordinary correlated random intercept
# shared by two responses.
#
# The known standard deviation arrives through `vreal(sd)` rather than
# `se(sd, sigma = FALSE)`. That is not a preference: the core gates
# `se()` on the family NAME ("se() is supported for gaussian and student
# families only"), so a custom family cannot opt in however faithfully
# it reads the term. `vreal()` is the channel a custom family IS given,
# and using it costs one word in the formula and one line here. The gap
# is recorded in dev/bcm-findings.md.
#
# There is no `sigma` parameter, because there is nothing to estimate: a
# free residual scale beside a known one is a flat direction.
bcm_gaussian_probit <- function() {
  frmtmb_family(
    "bcm_gaussian_probit",
    dpars = "mu",
    links = list(mu = bcm_probit()),
    type = "continuous",
    required_aterms = "vreal1",
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars[["mu"]], aterms[["vreal1"]], log = TRUE)
    },
    valid_y = function(y, aterms) {
      if (any(y < 0) || any(y > 1)) {
        stop("bcm_gaussian_probit(): the mean is Phi of a linear ",
             "predictor, so the response belongs on (0, 1); divide by ",
             "its scale first", call. = FALSE)
      }
      if (any(aterms[["vreal1"]] <= 0)) {
        stop("bcm_gaussian_probit(): vreal(sd) carries the known ",
             "measurement standard deviation, which must be positive",
             call. = FALSE)
      }
    },
    init_dpars = list(
      mu = function(y, aterms) min(max(mean(y), 0.02), 0.98)),
    post = list(
      mean_fn = function(dpars, aterms) dpars[["mu"]],
      var_fn = function(dpars, aterms) aterms[["vreal1"]]^2),
    sim = function(dpars, aterms, n) {
      stats::rnorm(n, dpars[["mu"]], aterms[["vreal1"]])
    })
}
