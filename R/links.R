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
         paste(names(frmtmb_links), collapse = ", "), call. = FALSE)
  }
  lk
}
