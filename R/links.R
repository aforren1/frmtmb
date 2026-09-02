# AD-safe link registry. stats::make.link clamps at C level in ways the AD
# tape cannot see, so link functions are defined explicitly over plain
# arithmetic that RTMB overloads.
#
# `logit_eta` and `log_eta` are the ROBUST fields, present only on the
# links that admit an exact one. An inverse link saturates: plogis(40) is
# exactly 1 in double precision and 1 - exp(-exp(4)) is too, so a density
# written over `1 - mu` reads log(0) = -Inf with an unusable gradient
# where the true log-density is a perfectly ordinary -40. The linear
# predictor never saturated, so these recover the log-odds (`logit_eta`)
# or the log mean (`log_eta`) from it directly, which is the quantity
# every such density actually needs. A link with no exact form leaves the
# field absent and the family keeps the plain round trip.
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

#' Look up a link by name in the AD-safe link registry. An already
#' resolved link list passes through unchanged, so callers can accept
#' either a name or a custom link. An unknown name errors and lists the
#' available links.
#'
#' @noRd
get_link <- function(name) {
  if (is.list(name)) return(name)
  lk <- frmtmb_links[[name]]
  if (is.null(lk)) {
    stop("Unknown link: '", name, "'. Available links: ",
         paste(names(frmtmb_links), collapse = ", "), call. = FALSE)
  }
  lk
}
