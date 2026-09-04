#' The Wiener first-passage time family
#'
#' A drift-diffusion model for a two-choice decision: a noisy evidence
#' accumulator starts between two boundaries and the response time is
#' the first time it touches one of them. The family models the response
#' time; which boundary was touched is data, supplied through the
#' `vint()` addition term.
#'
#' The parameterization is brms's `wiener` family, name for name:
#'
#' \describe{
#'   \item{`mu`}{Drift rate, the mean rate of evidence accumulation
#'     (brms and the literature also call it `v`). Identity link, so it
#'     is signed: positive drift favors the upper boundary.}
#'   \item{`bs`}{Boundary separation, the distance between the two
#'     boundaries (`a`). Log link.}
#'   \item{`ndt`}{Non-decision time, the part of the response time spent
#'     encoding and moving rather than deciding (`t0` or `tau`). Bounded
#'     link; see Non-decision time below.}
#'   \item{`bias`}{Relative start point in (0, 1), the fraction of the
#'     boundary separation the accumulator starts at (`w`). Logit link.
#'     0.5 is unbiased.}
#' }
#'
#' @section The decision indicator:
#' brms spells the boundary a trial ended at as `y | dec(decision)`.
#' frmtmb's formula grammar has no `dec()`, and its list of addition
#' terms is closed, so this family reads the indicator from `vint()`
#' instead:
#'
#' ```
#' frm(bf(rt | vint(upper) ~ condition), family = wiener(), data = dat)
#' ```
#'
#' `upper` must be 1 for a response at the upper boundary and 0 for the
#' lower one. Unlike brms's `dec()`, it will not accept a factor or the
#' strings `"upper"` and `"lower"`; code it yourself. Omitting `vint()`
#' altogether is refused with a message that says this, because the
#' failure is otherwise silent.
#'
#' @section Non-decision time:
#' The density is zero for a response time at or below `ndt`, so the
#' likelihood has a hard edge at `ndt = min(rt)` and an ordinary log
#' link would let the optimizer walk straight over it. The `ndt` link is
#' a logit scaled onto `(0, max_ndt)` instead, which makes the
#' constraint structural rather than a thing the optimizer has to
#' discover.
#'
#' `max_ndt` defaults to the smallest response time in the data, found
#' when the model frame is assembled. Give it explicitly to pin the
#' bound, which is worth doing when you will `predict()` on new data
#' whose minimum differs from the training minimum.
#'
#' Past a linear predictor of about 37 the logit saturates in double
#' precision and `ndt` rounds to `max_ndt` exactly. Nothing guards that,
#' and nothing needs to: the density falls off a cliff as the decision
#' time goes to zero, so the log likelihood is already unreachable long
#' before the link runs out of digits.
#'
#' @section Accuracy:
#' The density is the Navarro and Fuss (2009) pair of series, both
#' evaluated at a fixed truncation and combined with a smooth weight,
#' because an automatic-differentiation tape cannot choose between them
#' on a parameter. It agrees with `RWiener::dwiener()` to better than
#' 1e-12 relative on the log scale over normalized times from 1e-3 to
#' 50. See `vignette("ddm")`.
#'
#' @param max_ndt Upper bound for the non-decision time, in the units of
#'   the response. `NULL`, the default, takes it from the data.
#' @param link Link for the drift rate. Identity by default, and there
#'   is rarely a reason to change it: the drift rate is signed.
#'
#' @return A `frmtmb_family`.
#'
#' @references
#' Navarro, D. J. and Fuss, I. G. (2009). Fast and accurate calculations
#' for first-passage times in Wiener diffusion models. *Journal of
#' Mathematical Psychology*, 53(4), 222-230.
#'
#' @examples
#' set.seed(1)
#' dat <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
#' fit <- frm(bf(rt | vint(upper) ~ 1, bias = 0.5),
#'            family = wiener(), data = dat)
#' fixef(fit)
#'
#' @export
wiener <- function(max_ndt = NULL, link = "identity") {
  if (!is.null(max_ndt)) {
    if (!is.numeric(max_ndt) || length(max_ndt) != 1L ||
        !is.finite(max_ndt) || max_ndt <= 0) {
      stop("wiener(): `max_ndt` must be one positive finite number, ",
           "or NULL to take it from the data.", call. = FALSE)
    }
  }
  # The bound lives in an environment the link closures read at call
  # time, so that a default bound can be filled in from the response
  # once frmtmb has the data. valid_y() is the only slot that sees the
  # response before a link is used.
  st <- new.env(parent = emptyenv())
  st$ub <- if (is.null(max_ndt)) NA_real_ else max_ndt
  st$fixed <- !is.null(max_ndt)

  bound <- function() {
    if (is.na(st$ub)) {
      stop("wiener(): the non-decision time bound is not set yet. This ",
           "happens when a wiener() family object is used outside frm(), ",
           "for example to inspect its links before a fit. Pass ",
           "`max_ndt` to wiener() to set the bound up front.",
           call. = FALSE)
    }
    st$ub
  }

  ndt_link <- list(
    name = "scaled_logit",
    linkfun = function(mu) { U <- bound(); log(mu / (U - mu)) },
    linkinv = function(eta) bound() / (1 + exp(-eta)),
    mu_eta = function(eta) {
      p <- 1 / (1 + exp(-eta))
      bound() * p * (1 - p)
    })

  frmtmb::custom_family(
    "wiener",
    dpars = c("mu", "bs", "ndt", "bias"),
    links = list(mu = link, bs = "log", ndt = ndt_link, bias = "logit"),
    lpdf = function(y, dpars, aterms) {
      up <- aterms[["vint1"]]
      ddm_lpdf_both(y - dpars[["ndt"]],
                        dpars[["mu"]], dpars[["bs"]], dpars[["bias"]], up)
    },
    valid_y = function(y, aterms) {
      ddm_check_response(y, aterms, st)
    },
    init_dpars = list(
      # A drift of zero is the honest starting guess: the sign of the
      # drift is what the data are there to tell us, and a start with
      # the wrong sign costs more than a start at the middle.
      mu = function(y, aterms) 0,
      bs = function(y, aterms) 1.5,
      ndt = function(y, aterms) 0.5 * min(y),
      bias = function(y, aterms) 0.5),
    type = "continuous",
    post = list(mean_fn = function(dpars, aterms) {
      ddm_mean_rt(dpars, aterms)
    }),
    sim = function(dpars, aterms, n) {
      ddm_sim_rt(dpars, aterms, n)
    },
    sim_refusal = NULL)
}

#' Response and decision-indicator validation.
#'
#' Also the hook that fills in the default non-decision-time bound:
#' `valid_y()` is called once when the model frame is assembled, and
#' before any link function runs, so it is the one place a family can
#' learn something about the response and act on it.
#'
#' @noRd
ddm_check_response <- function(y, aterms, st) {
  if (any(!is.finite(y)) || any(y <= 0)) {
    stop("wiener: the response must be a strictly positive, finite ",
         "response time.", call. = FALSE)
  }
  up <- aterms[["vint1"]]
  if (is.null(up)) {
    # frmtmb has no way for a family to declare that it requires an
    # addition term, and a missing vint() would otherwise reach the
    # density as NULL, where the arithmetic silently collapses to a
    # zero-length log-likelihood and the fit "succeeds".
    stop("wiener: the decision indicator is missing. Which boundary a ",
         "trial ended at is data, and it reaches the family through ",
         "vint():\n",
         "    frm(bf(rt | vint(upper) ~ x), family = wiener(), ...)\n",
         "where `upper` is 1 for the upper boundary and 0 for the ",
         "lower one. brms spells this `rt | dec(decision)`; frmtmb's ",
         "addition terms are a closed set and dec() is not one of them.",
         call. = FALSE)
  }
  if (any(!is.finite(up)) || any(up != 0 & up != 1)) {
    stop("wiener: vint() decision indicator must be 0 (lower boundary) ",
         "or 1 (upper boundary). brms's dec() accepts a factor or the ",
         "strings \"upper\"/\"lower\"; this one does not, so recode it ",
         "with as.integer(decision == \"upper\").", call. = FALSE)
  }
  if (is.na(st$ub)) {
    st$ub <- min(y)
  } else if (st$fixed && st$ub > min(y)) {
    # min(y) itself is allowed, and is the default: the scaled logit
    # never reaches its own bound at a finite linear predictor, so
    # ndt < min(rt) stays strict. Anything above min(y) does admit
    # parameter values with no likelihood.
    stop("wiener: max_ndt = ", format(st$ub), " is above the smallest ",
         "response time (", format(min(y)), "). The density is zero at ",
         "and below the non-decision time, so a bound above min(rt) ",
         "admits parameter values with no likelihood.", call. = FALSE)
  }
  invisible(NULL)
}
