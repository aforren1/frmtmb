#' Link functions
#'
#' The links frmtmb accepts, what each one maps, and where each one may
#' be used. Every family constructor takes a `link` argument for the
#' mean and a `link_<dpar>` argument for each of its other
#' distributional parameters, so `student(link_sigma = "softplus")`
#' puts sigma on a softplus. The roster follows brms 2.23.0: a model
#' ported from brms is fitted through the same inverse link.
#'
#' A link is defined here over plain arithmetic that RTMB overloads,
#' not through [stats::make.link()], because `make.link()` clamps at C
#' level in ways the automatic differentiation tape cannot see.
#'
#' @section The roster:
#'
#' `g` is the link and `g^-1` its inverse. The inverse is the one that
#' does the work: it maps a linear predictor to the parameter, so its
#' RANGE is the parameter support the link fits. `eta` is the linear
#' predictor and `mu` the parameter.
#'
#' | Link | `g(mu)` | `g^-1(eta)` | Range | Robust field |
#' |---|---|---|---|---|
#' | `identity` | `mu` | `eta` | the real line | |
#' | `log` | `log(mu)` | `exp(eta)` | positive | `log_eta` |
#' | `softplus` | `log(expm1(mu))` | `log1p(exp(eta))` | positive | `log_eta` |
#' | `squareplus` | `(mu^2 - 1) / mu` | `(eta + sqrt(eta^2 + 4)) / 2` | positive | `log_eta` |
#' | `inverse` | `1 / mu` | `1 / eta` | positive | |
#' | `1/mu^2` | `1 / mu^2` | `1 / sqrt(eta)` | positive | |
#' | `sqrt` | `sqrt(mu)` | `eta^2` | non-negative | `log_eta` |
#' | `logit` | `log(mu / (1 - mu))` | `1 / (1 + exp(-eta))` | the unit interval | `logit_eta` |
#' | `probit` | `qnorm(mu)` | `pnorm(eta)` | the unit interval | `logit_eta` |
#' | `probit_approx` | `qnorm(mu)` | `plogis(0.07056 eta^3 + 1.5976 eta)` | the unit interval | `logit_eta` |
#' | `cloglog` | `log(-log(1 - mu))` | `1 - exp(-exp(eta))` | the unit interval | `logit_eta` |
#' | `cauchit` | `tan(pi (mu - 0.5))` | `0.5 + atan(eta) / pi` | the unit interval | |
#' | `softit` | `log(expm1(mu / (1 - mu)))` | `s / (1 + s)` for `s = log1p(exp(eta))` | the unit interval | `logit_eta` |
#' | `logm1` | `log(mu - 1)` | `1 + exp(eta)` | above one | |
#' | `log1p` | `log1p(mu)` | `expm1(eta)` | above minus one | |
#' | `power12` | `log((mu - 1) / (2 - mu))` | `1 + 1 / (1 + exp(-eta))` | between one and two | |
#' | `tan_half` | `tan(mu / 2)` | `2 atan(eta)` | the circle | |
#'
#' `probit_approx` is the one pair that is not an exact inverse, in
#' frmtmb and in brms alike. brms generates Stan code that uses
#' `Phi_approx()`, the logistic of a cubic, while its own R-side
#' `inv_link()` answers `pnorm()`. The Stan form is the one a ported
#' model was fitted with, so it is the form used here.
#'
#' `sqrt` is not injective on the whole line. brms allows it for a
#' count mean anyway, and so does frmtmb.
#'
#' @section Links for the mean:
#'
#' Any link in the roster is accepted for `mu`. The table below is what
#' brms 2.23.0 accepts, so it says which pairings PORT. frmtmb does not
#' refuse the others, because an extension family is free to mean
#' something else by its own `mu`. The first link listed is the
#' default.
#'
#' | Family | Links for `mu` |
#' |---|---|
#' | `gaussian`, `student`, `skew_normal` | `identity`, `log`, `inverse`, `softplus`, `squareplus`, and every unit-interval link |
#' | `lognormal`, `shifted_lognormal`, `hurdle_lognormal` | `identity`, `inverse` |
#' | `Gamma`, `weibull`, `exponential`, `hurdle_gamma` | `log`, `identity`, `inverse`, `softplus`, `squareplus` |
#' | `asym_laplace`, `exgaussian`, `zero_inflated_asym_laplace` | `identity`, `log`, `inverse`, `softplus`, `squareplus` |
#' | `poisson`, `negbinomial`, `geometric`, `compois`, and the zero-inflated and hurdle counts | `log`, `identity`, `sqrt`, `softplus`, `squareplus` |
#' | `binomial`, `bernoulli`, `Beta`, `zero_inflated_binomial`, `zero_inflated_beta` | `logit`, `probit`, `probit_approx`, `cloglog`, `cauchit`, `softit`, `identity`, `log` |
#' | `beta_binomial` | the same list WITHOUT `log`: brms's own table for this one family stops at `identity` |
#' | `inverse.gaussian` | `1/mu^2`, `inverse`, `identity`, `log`, `softplus`, `squareplus` |
#' | `cox` | `log`, `identity`, `softplus`, `squareplus` |
#' | `von_mises` | `tan_half`, `identity` |
#' | `cumulative` | `logit`, `probit`, `probit_approx`, `cloglog`, `cauchit`, `softit` |
#' | `sratio`, `cratio` | `logit`, `probit`, `probit_approx`, `cloglog`, `cauchit` |
#' | `acat` | `logit` ONLY. brms takes the same six as `cumulative`; this is the one place frmtmb departs, and the reason is below |
#' | `categorical`, `multinomial` | `logit` |
#'
#' An ordinal family is the exception that IS enforced. Its `link`
#' names the cumulative distribution function the thresholds are read
#' through, not a link on a mean, so only a link whose inverse maps
#' onto the unit interval can serve. The rest are refused by name.
#'
#' `acat()` is refused for a different reason, and it is the one place
#' in this table where frmtmb takes less than brms. Probit, cloglog,
#' cauchit and softit all map onto the unit interval and brms accepts
#' every one of them for `acat`, so the refusal is not about the link.
#' It is about the density. `brms:::inv_link_acat()` branches: on the
#' logit a category probability is `c(1, cumprod(exp(x)))` normalized,
#' which is the log-linear form frmtmb implements, and off the logit it
#' is a product of distribution functions times a reversed product of
#' survivals, normalized. The second form agrees with the first when
#' the distribution function is logistic, so it generalizes the same
#' model rather than replacing it, but it is a second expression that
#' has to be written and taped. Substituting a distribution function
#' into the log-linear form does not reach it, which is why the other
#' three ordinal families could be routed through this registry and
#' `acat()` could not. `acat()` says all of this when it refuses.
#'
#' @section Links for the other distributional parameters:
#'
#' `link_<dpar>` takes the values brms allows for that parameter, and
#' the set is fixed by what the parameter must stay inside. The first
#' link listed is the default.
#'
#' | Parameter | Argument | Links |
#' |---|---|---|
#' | `sigma`, `shape`, `phi`, `kappa`, `beta`, `ndt`, and the `nu` of `compois` | `link_sigma` and so on | `log`, `identity`, `softplus`, `squareplus` |
#' | `nu` of `student` | `link_nu` | `logm1`, `identity` |
#' | `zi`, `hu`, `quantile` | `link_zi`, `link_hu`, `link_quantile` | `logit`, `identity` |
#' | `alpha` of `skew_normal` | `link_alpha` | `identity`, `log`, `softplus`, `squareplus` |
#'
#' `identity` is in every one of those sets because brms puts it there.
#' It lets a parameter that must stay positive go negative, which the
#' density then reports as a `NaN` rather than as a bad link. It is
#' offered so that a brms model ports, not because it is a good choice.
#'
#' The `power` of `tweedie()` takes no argument. It is confined to the
#' open interval from one to two, `power12` is the only link in the
#' roster that maps onto it, and brms has no tweedie to port from.
#'
#' A prior class is a distributional parameter name, so changing a
#' parameter's link changes where its prior sits. `set_prior(class =
#' "sigma")` is a density on sigma itself, and the placement `frm()`
#' gives it is read from sigma's link, so changing `link_sigma` changes
#' the Jacobian that prior carries. See [set_prior()].
#'
#' @section What a robust field buys:
#'
#' An inverse link SATURATES. `plogis(40)` is exactly one in double
#' precision, and so is `1 - exp(-exp(4))`. A density written over
#' `1 - mu` then reads `log(0)`, which is `-Inf` with an unusable
#' gradient, where the true log density is an ordinary `-40`. The
#' linear predictor never saturated, so a link that carries a robust
#' field recovers the quantity the density actually needs from the
#' linear predictor directly: `logit_eta` gives the log odds and
#' `log_eta` the log mean.
#'
#' The field is present only where the plain round trip measurably
#' saturates AND an exact form exists. `cauchit` has neither the
#' problem nor a cure: its tails are polynomial, so the plain round
#' trip keeps eight digits out to a linear predictor of 3.2e9.
#' `inverse` and `1/mu^2` cannot overflow or underflow at any
#' representable positive linear predictor.
#'
#' Nothing has to be done to use this. A family whose link carries the
#' field is put on the robust form of its density automatically:
#' `binomial()` on `dbinom_robust()`, `negbinomial()` on
#' `dnbinom_robust()`.
#'
#' @section Custom links:
#'
#' Wherever a link is taken, a list is taken instead of a name. It must
#' carry `name`, `linkfun`, `linkinv` and `mu_eta`, the derivative of
#' `linkinv`, which `predict(se.fit = TRUE)` and every delta-method
#' interval read. `logit_eta` and `log_eta` are optional; supply one
#' only if it is exact.
#'
#' @return A link is not a free-standing function. It is named to a
#'   family constructor, which resolves it against the registry and
#'   carries it on the fitted model, so the value a link contributes is
#'   the scale each distributional parameter is estimated on. Read it
#'   back with `frm_family()`, and read a custom link's own fields with
#'   the list that was supplied. This page documents the roster.
#' @seealso [frmtmb-families] for the constructors that take these,
#'   [set_prior()] for what a link does to a prior.
#' @name frmtmb-links
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(80))
#' d$y <- rbinom(80, 1, pnorm(0.4 + 0.8 * d$x))
#'
#' # the mean on a probit rather than a logit
#' fixef(frm(bf(y ~ x), family = bernoulli(link = "probit"), data = d))$mu
#'
#' # a link on a parameter that is not the mean
#' d$z <- rnorm(80, 1 + d$x, exp(0.2 + 0.3 * d$x))
#' frm(bf(z ~ x, sigma ~ x), family = student(link_sigma = "log"), data = d)
NULL

# AD-safe link registry. stats::make.link clamps at C level in ways the AD
# tape cannot see, so link functions are defined explicitly over plain
# arithmetic that RTMB overloads.
#
# The roster follows brms, so that a model ported from it is fitted
# through the same inverse link. `pnorm` and `qnorm` are written
# RTMB-qualified because the `stats` versions have no advector method;
# everything else here is a base generic RTMB overloads.
#
# `logit_eta` and `log_eta` are the ROBUST fields, present only on the
# links whose plain round trip measurably saturates AND that admit an
# exact form. An inverse link saturates: plogis(40) is exactly 1 in
# double precision and 1 - exp(-exp(4)) is too, so a density written
# over `1 - mu` reads log(0) = -Inf with an unusable gradient where the
# true log-density is a perfectly ordinary -40. The linear predictor
# never saturated, so these recover the log-odds (`logit_eta`) or the
# log mean (`log_eta`) from it directly, which is the quantity every
# such density actually needs. A link with no exact form, and a link
# whose tails are polynomial rather than exponential and so never
# saturate anywhere an optimizer can reach, leaves the field absent and
# the family keeps the plain round trip. Each entry below states which
# it is and at what linear predictor the measurement put the boundary.
frmtmb_links <- list(
  identity = list(
    name    = "identity",
    linkfun = function(mu) mu,
    linkinv = function(eta) eta,
    mu_eta  = function(eta) rep(1, length(eta))
  ),
  log = list(
    name    = "log",
    linkfun = function(mu) log(mu),
    linkinv = function(eta) exp(eta),
    mu_eta  = function(eta) exp(eta),
    log_eta = function(eta) eta
  ),
  logit = list(
    name    = "logit",
    linkfun = function(mu) log(mu / (1 - mu)),
    linkinv = function(eta) 1 / (1 + exp(-eta)),
    mu_eta  = function(eta) {
      p <- 1 / (1 + exp(-eta))
      p * (1 - p)
    },
    logit_eta = function(eta) eta
  ),
  cloglog = list(
    name    = "cloglog",
    linkfun = function(mu) log(-log(1 - mu)),
    linkinv = function(eta) 1 - exp(-exp(eta)),
    mu_eta  = function(eta) exp(eta - exp(eta)),
    # log(mu) - log(1 - mu) with log(1 - mu) = -exp(eta) exactly and
    # log(mu) = log(1 - exp(-exp(eta))) through expm1, which keeps the
    # small-mu end. cloglog saturates far earlier than the logit - at
    # eta = 4 the round trip already gives mu = 1 - so this matters at
    # single-digit linear predictors.
    logit_eta = function(eta) {
      t <- exp(eta)
      log(-expm1(-t)) + t
    }
  ),
  # The standard normal CDF: probit regression, and signal detection
  # theory whenever it is written as a model rather than as d-prime.
  probit = list(
    name    = "probit",
    linkfun = function(mu) RTMB::qnorm(mu),
    linkinv = function(eta) RTMB::pnorm(eta),
    # the normal density, written out rather than taken from dnorm:
    # RTMB::pnorm's own AD derivative is this to the last bit
    mu_eta  = function(eta) exp(-0.5 * eta^2) / sqrt(2 * pi),
    # 1 - pnorm(eta) rounds to zero once pnorm(eta) rounds to one, and
    # the plain log(1 - mu) is past 1e-8 relative accuracy by eta = 6.4.
    # A difference of two pnorm logs cannot cancel, and holds until
    # pnorm itself underflows at |eta| = 38.2 - which is also where the
    # plain log(mu) gives out, so the other tail loses nothing.
    logit_eta = function(eta) {
      log(RTMB::pnorm(eta)) - log(RTMB::pnorm(-eta))
    }
  ),
  # brms's probit_approx. Its Stan program uses Phi_approx(), which is
  # inv_logit(0.07056 x^3 + 1.5976 x), while brms's own R-side
  # inv_link() answers pnorm() for this link; the two disagree, and the
  # Stan form is the one a ported model was fitted with. linkfun is
  # qnorm, brms's choice, so the pair is not an exact inverse in either
  # package.
  probit_approx = list(
    name    = "probit_approx",
    linkfun = function(mu) RTMB::qnorm(mu),
    linkinv = function(eta) {
      1 / (1 + exp(-(0.07056 * eta^3 + 1.5976 * eta)))
    },
    mu_eta  = function(eta) {
      p <- 1 / (1 + exp(-(0.07056 * eta^3 + 1.5976 * eta)))
      p * (1 - p) * (0.21168 * eta^2 + 1.5976)
    },
    # the inverse link is a logistic OF the cubic, so the log-odds is
    # the cubic exactly and never saturates at any finite eta. The plain
    # round trip is past 1e-8 by eta = 5.9.
    logit_eta = function(eta) 0.07056 * eta^3 + 1.5976 * eta
  ),
  # The standard Cauchy CDF. No logit_eta: the Cauchy tail is
  # polynomial, not exponential - 1 - mu is 1 / (pi eta) to leading
  # order - so the plain round trip keeps 1e-8 relative accuracy out to
  # eta = 3.2e9 and there is nothing for a robust field to rescue. Nor
  # would one be writable: the exact upper tail atan(1 / eta) / pi holds
  # only above zero, and below it the same expression is off by pi.
  cauchit = list(
    name    = "cauchit",
    linkfun = function(mu) tan(pi * (mu - 0.5)),
    linkinv = function(eta) 0.5 + atan(eta) / pi,
    mu_eta  = function(eta) 1 / (pi * (1 + eta^2))
  ),
  # brms's softit: a softplus squashed onto (0, 1).
  softit = list(
    name    = "softit",
    # log(expm1(r)) as r + log(-expm1(-r)), the same number without
    # expm1 overflowing at r = 710
    linkfun = function(mu) {
      r <- mu / (1 - mu)
      r + log(-expm1(-r))
    },
    linkinv = function(eta) {
      y <- RTMB::logspace_add(0 * eta, eta)
      y / (1 + y)
    },
    mu_eta  = function(eta) {
      y <- RTMB::logspace_add(0 * eta, eta)
      1 / ((1 + exp(-eta)) * (1 + y)^2)
    },
    # mu / (1 - mu) is the softplus exactly, so the log-odds is its log.
    # The plain round trip loses the upper tail at eta = 3.2e10, where
    # y / (1 + y) rounds to one; both paths end together at eta = -745,
    # where the softplus itself underflows.
    logit_eta = function(eta) log(RTMB::logspace_add(0 * eta, eta))
  ),
  inverse = list(
    name    = "inverse",
    linkfun = function(mu) 1 / mu,
    linkinv = function(eta) 1 / eta,
    mu_eta  = function(eta) -1 / eta^2
  ),
  # log(x - 1): keeps student-t df above 1 (brms convention for nu)
  logm1 = list(
    name    = "logm1",
    linkfun = function(mu) log(mu - 1),
    linkinv = function(eta) 1 + exp(eta),
    mu_eta  = function(eta) exp(eta)
  ),
  # log(x + 1): brms offers it on one dpar only, the xi of
  # gen_extreme_value, whose support is (-1, Inf). Neither robust field
  # is meaningful for a dpar that may be negative.
  log1p = list(
    name    = "log1p",
    linkfun = function(mu) log1p(mu),
    linkinv = function(eta) expm1(eta),
    mu_eta  = function(eta) exp(eta)
  ),
  # tan of the half angle: maps the whole line onto the circle
  # (-pi, pi), which is the support of a von Mises mean direction
  # (brms's tan_half link)
  tan_half = list(
    name    = "tan_half",
    linkfun = function(mu) tan(mu / 2),
    linkinv = function(eta) 2 * atan(eta),
    mu_eta  = function(eta) 2 / (1 + eta^2)
  ),
  # logit onto (1, 2): the tweedie power parameter's valid range
  power12 = list(
    name    = "power12",
    linkfun = function(mu) log((mu - 1) / (2 - mu)),
    linkinv = function(eta) 1 + 1 / (1 + exp(-eta)),
    mu_eta  = function(eta) {
      p <- 1 / (1 + exp(-eta))
      p * (1 - p)
    }
  ),
  # brms's softplus, log1p(exp(eta)), for a strictly positive dpar.
  # Through logspace_add rather than log(1 + exp(eta)), which overflows
  # at eta = 710 where the answer is an ordinary 710.
  softplus = list(
    name    = "softplus",
    # log(expm1(mu)) without the overflow, as in softit above
    linkfun = function(mu) mu + log(-expm1(-mu)),
    linkinv = function(eta) RTMB::logspace_add(0 * eta, eta),
    mu_eta  = function(eta) 1 / (1 + exp(-eta)),
    # Exact wherever the softplus itself is, which is down to
    # eta = -745. What this buys is not accuracy in log(mu) - the round
    # trip has that already - but the density BRANCH. Without it
    # negbinomial forms `mu + mu^2 / shape`, and at mu = 9e-14 that sum
    # has lost the excess over mu to rounding, so the log-density
    # carries noise that a central difference sees and the tape does
    # not. log_eta is what puts the family on dnbinom_robust instead.
    log_eta = function(eta) log(RTMB::logspace_add(0 * eta, eta))
  ),
  # brms's squareplus, for a strictly positive dpar.
  squareplus = list(
    name    = "squareplus",
    linkfun = function(mu) (mu^2 - 1) / mu,
    linkinv = function(eta) (eta + sqrt(eta^2 + 4)) / 2,
    mu_eta  = function(eta) (1 + eta / sqrt(eta^2 + 4)) / 2,
    # (eta + sqrt(eta^2 + 4)) / 2 IS exp(asinh(eta / 2)), so the log
    # mean is asinh(eta / 2) exactly - and asinh has none of the
    # cancellation the sum has. Below eta = -1e5 the plain form has lost
    # 1e-8 of the tail, and it eventually returns zero, where asinh goes
    # on answering -log(-eta).
    log_eta = function(eta) asinh(eta / 2)
  ),
  # brms's sqrt link for a count mean. Not injective on the line; brms
  # allows it anyway and so does this.
  sqrt = list(
    name    = "sqrt",
    linkfun = function(mu) sqrt(mu),
    linkinv = function(eta) eta^2,
    mu_eta  = function(eta) 2 * eta,
    # log(eta^2) rather than 2 * log(abs(eta)), which is the same number
    # without a kink on the tape at eta = 0, and exact until eta^2
    # underflows at |eta| = 1e-160. Present for the same reason as
    # softplus's: it selects the robust density branch, which is what
    # keeps a count family honest at a mean near zero.
    log_eta = function(eta) log(eta^2)
  ),
  # The inverse Gaussian canonical link, in brms's spelling. Unrelated
  # to power12 above, which is the tweedie logit onto (1, 2). No
  # log_eta: 1 / sqrt(eta) can neither underflow nor overflow at any
  # representable positive eta, so the plain log(mu) never loses a digit
  # and -0.5 * log(eta) would buy nothing.
  `1/mu^2` = list(
    name    = "1/mu^2",
    linkfun = function(mu) 1 / mu^2,
    linkinv = function(eta) 1 / sqrt(eta),
    mu_eta  = function(eta) -0.5 * eta^(-1.5)
  )
)

#' `log(p)` from a log-odds, as `-log(1 + exp(-x))`. RTMB's
#' `logspace_add()` is exact at both ends and differentiable through
#' them, which `log(plogis(x))` is not: it loses the whole upper tail to
#' rounding.
#'
#' @noRd
log_inv_logit <- function(x) -RTMB::logspace_add(0 * x, -x)

#' `log(1 - p)` from a log-odds. The mirror of [log_inv_logit()], and the
#' term every saturating binomial-style density is actually missing.
#'
#' @noRd
log1m_inv_logit <- function(x) -RTMB::logspace_add(0 * x, x)

#' The fields every link object must carry. `linkfun` and `linkinv` move
#' between the scales; `mu_eta` is the derivative the delta method needs.
#' The robust fields (`logit_eta`, `log_eta`) are optional by design.
#'
#' @noRd
link_required_fields <- c("name", "linkfun", "linkinv", "mu_eta")

#' Look up a link by name in the AD-safe link registry. An already
#' resolved link list passes through, after the four required fields are
#' checked: a custom link used to be accepted untouched, and one missing
#' `mu_eta` fit, summarized and predicted happily before failing inside
#' `predict(se.fit = TRUE)`, a call site with nothing to say about the
#' family that caused it. `dpar` names that family slot when there is
#' one. An unknown name errors and lists the available links.
#'
#' @noRd
get_link <- function(name, dpar = NULL) {
  if (is.list(name)) {
    where <- if (is.null(dpar)) "" else paste0(" of dpar '", dpar, "'")
    absent <- setdiff(link_required_fields, names(name))
    if (length(absent)) {
      stop("The custom link", where, " has no ",
           paste0("`", absent, "`", collapse = ", "),
           ". A link object needs name, linkfun, linkinv and mu_eta ",
           "(the derivative of linkinv, which predict(se.fit = TRUE) ",
           "and every delta-method interval read)", call. = FALSE)
    }
    bad <- Filter(function(f) !is.function(name[[f]]),
                  c("linkfun", "linkinv", "mu_eta"))
    if (length(bad)) {
      stop("The custom link", where, " has a non-function ",
           paste0("`", bad, "`", collapse = ", "),
           "; each must be a function of one vector", call. = FALSE)
    }
    # [[ ]]: `$` on a link list is how a partial match would silently
    # answer for a field that is not there
    nm_field <- name[["name"]]
    if (!is.character(nm_field) || length(nm_field) != 1L ||
        is.na(nm_field)) {
      stop("The custom link", where, " must name itself with a single ",
           "string in `name`; it labels the link in summary() and in ",
           "every method that reports the scale", call. = FALSE)
    }
    return(name)
  }
  # `[[` on a list indexes RECURSIVELY when given a vector, so
  # frmtmb_links[[c("log", "name")]] is frmtmb_links$log$name, the string
  # "log", which is not a link at all and was handed back as one. An
  # integer index picks a link by position, so link = 1L silently became
  # the identity link. Neither reaches the "Unknown link" branch.
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    stop("A link must be named by a single string, e.g. link = \"logit\", ",
         "not ", arg_desc(name), ". Available links: ",
         paste(names(frmtmb_links), collapse = ", "), call. = FALSE)
  }
  lk <- frmtmb_links[[name]]
  if (is.null(lk)) {
    stop("Unknown link: '", name, "'. Available links: ",
         paste(names(frmtmb_links), collapse = ", "),
         ". See ?`frmtmb-links` for what each one maps and which ",
         "families take it", call. = FALSE)
  }
  lk
}

# The link sets a distributional parameter other than the mean may
# take, measured from brms:::links_dpars() in brms 2.23.0. The FIRST
# entry is brms's default and frmtmb's.
#
# Keyed by the parameter's SUPPORT and passed in from the call site,
# not looked up by the parameter's name. A name lookup gets `nu` wrong:
# brms answers `logm1` for it, which is student's degrees of freedom,
# but frmtmb's compois also calls its dispersion `nu` and that is an
# ordinary positive scale. The constructor knows which one it has.
dpar_links_positive <- c("log", "identity", "softplus", "squareplus")
dpar_links_unit     <- c("logit", "identity")
dpar_links_above1   <- c("logm1", "identity")
# alpha of skew_normal: the same four as a positive parameter, but the
# skewness is signed, so identity is the default rather than a way to
# break the density
dpar_links_signed   <- c("identity", "log", "softplus", "squareplus")

#' Resolve a `link_<dpar>` constructor argument to a link object.
#'
#' `choices` is what the parameter's support admits, so an out-of-range
#' link is refused HERE, with the parameter named, rather than at the
#' first `NaN` in the density where nothing says which link caused it.
#' A list still passes through to [get_link()], so an extension family
#' can supply a custom link for a dpar exactly as it can for the mean.
#'
#' @noRd
dpar_link <- function(value, dpar, family, choices) {
  if (is.list(value)) return(get_link(value, dpar = dpar))
  allowed <- paste0("\"", choices, "\"", collapse = ", ")
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop(family, "(link_", dpar, " =) takes a single link name, not ",
         arg_desc(value), ". Allowed: ", allowed,
         ". See ?`frmtmb-links`", call. = FALSE)
  }
  if (!value %in% choices) {
    stop(family, "(link_", dpar, " =): '", value, "' is not a link `",
         dpar, "` can take. Allowed: ", allowed,
         ", the set brms allows, fixed by the range `", dpar,
         "` has to stay inside. See ?`frmtmb-links`", call. = FALSE)
  }
  get_link(value, dpar = dpar)
}
