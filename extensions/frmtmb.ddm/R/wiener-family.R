#' The Wiener first-passage time family
#'
#' A drift-diffusion model for a two-choice decision: a noisy evidence
#' accumulator starts between two boundaries and the response time is
#' the first time it touches one of them. The family models the response
#' time; which boundary was touched is data, supplied through the
#' `dec()` addition term.
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
#' `variability` adds Ratcliff's three across-trial variability
#' parameters to that set; see Across-trial variability below.
#'
#' @section The decision indicator:
#' The boundary a trial ended at is data, and it reaches the density as
#' an addition term:
#'
#' ```
#' frm(bf(rt | dec(response) ~ condition), family = wiener(), data = dat)
#' ```
#'
#' `dec()` is spelled as brms spells it and takes what brms takes: a
#' factor or character vector whose SECOND level is the upper boundary
#' (so `"lower"`/`"upper"` and `c(FALSE, TRUE)` both work as they read),
#' or a numeric 0/1 column. The package contributes the term to frmtmb's
#' addition-term registry when it loads.
#'
#' `vint()` also still works, and carries the indicator as a plain 0/1
#' integer column:
#'
#' ```
#' frm(bf(rt | vint(upper) ~ condition), family = wiener(), data = dat)
#' ```
#'
#' Supplying neither is refused with a message that says so, because the
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
#' A `max_ndt` above the smallest response time is refused, because for
#' this family alone it admits parameter values at which some observed
#' row has no likelihood. Inside a [frmtmb::mixture()] that is exactly
#' what the other component is for, so `allow_unreachable = TRUE` lifts
#' the refusal; see Mixtures.
#'
#' Past a linear predictor of about 37 the logit saturates in double
#' precision and `ndt` rounds to `max_ndt` exactly. Nothing guards that,
#' and nothing needs to: the density falls off a cliff as the decision
#' time goes to zero, so the log likelihood is already unreachable long
#' before the link runs out of digits.
#'
#' @section Across-trial variability:
#' Ratcliff's full diffusion model draws three of the four parameters
#' afresh on every trial. `variability` names which of those to
#' estimate, and each one it names becomes an ordinary distributional
#' parameter that takes its own formula:
#'
#' ```
#' frm(bf(rt | dec(response) ~ coherence, bias = 0.5),
#'     family = wiener(variability = c("sv", "sz", "st")), data = dat)
#' ```
#'
#' \describe{
#'   \item{`sv`}{Standard deviation of a normal drift rate. Log link.}
#'   \item{`sz`}{Width of a uniform relative start point, centered on
#'     `bias`, on the same (0, 1) scale as `bias`. Logit link, so the
#'     width is below 1 and the start point stays inside the boundaries
#'     whenever `bias` is 0.5.}
#'   \item{`st`}{Width of a uniform non-decision time, centered on
#'     `ndt`, in the units of the response. Logit link scaled onto
#'     `(0, 2 * max_ndt)`.}
#' }
#'
#' The likelihood is the analytic Wiener density averaged over those
#' distributions, and the three are done three different ways because
#' they are three different integrals. The drift integral is Gaussian
#' against an exponential-quadratic and is evaluated in CLOSED FORM: it
#' is exact, it takes no nodes, and there is nothing to tune. The other
#' two are uniform and are evaluated by fixed-node Gauss-Legendre
#' quadrature, whose node counts are the `nodes` argument.
#'
#' The node positions and counts are decided when the family object is
#' built and are constants from then on, because a node count that moved
#' with a parameter would be a branch on a parameter and an
#' automatic-differentiation tape cannot record one. A parameter only
#' rescales the interval the fixed nodes are mapped onto.
#'
#' `fitted()` and `simulate()` follow the variability rather than
#' ignoring it. Both condition on the boundary a row ended at, and
#' conditioning reweights which per-trial parameters that row could have
#' had: the trials that reached a boundary are not a fair sample of the
#' drift rates that could have produced them. So the fitted mean is a
#' ratio of two quadratures and the simulator accepts a drawn drift rate
#' and start point with the boundary probability they imply. The plain
#' closed-form mean is not a usable approximation here: at an unbiased
#' start point it returns the same number for both boundaries, and the
#' full model does not.
#'
#' Two limits of the parameterization are worth knowing. The uniform
#' start point stays inside the boundaries by construction only when
#' `bias` is 0.5, which is the usual case and the one the logit link on
#' `sz` is scaled for; at a strongly biased start a wide `sz` can push
#' the range past a boundary, where the density is near zero and the
#' likelihood is a barrier rather than a cliff. And nothing holds
#' `ndt - st / 2` above zero, so a fit is free to report a range that
#' includes a negative non-decision time. Neither can be made structural
#' from outside frmtmb: both are joint constraints on two distributional
#' parameters, and a link is a property of one.
#'
#' @section Mixtures:
#' A contaminant component covers the trials the diffusion process
#' cannot produce, which is the standard treatment for fast guesses. The
#' Wiener component then wants a non-decision time that some observed
#' rows fall below, so pass the bound and lift the refusal:
#'
#' ```
#' frm(bf(rt | dec(response) ~ 1, bias1 = 0.5),
#'     family = mixture(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
#'                      lognormal()),
#'     data = dat)
#' ```
#'
#' A row below the non-decision time gets a log density of about
#' `-1 / delta` where `delta` is a billionth of the smallest response
#' time: it exponentiates to exactly zero, which is the right likelihood
#' for the component, and it differentiates to exactly zero, which a
#' true `-Inf` would not.
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
#' @param variability Which across-trial variability parameters to
#'   estimate: any of `"sv"` (drift rate), `"sz"` (start point) and
#'   `"st"` (non-decision time). The default estimates none, which is
#'   the plain Wiener model.
#' @param nodes Gauss-Legendre node counts for the two quadratures, as a
#'   named vector. Only the entries for the `variability` parameters in
#'   use are read. The defaults are measured rather than chosen: `sz`
#'   reaches 1e-10 in the log density by 7 nodes everywhere it was
#'   probed, and `st` needs more because a response time below
#'   `ndt + st / 2` cuts the range and the integrand turns on sharply at
#'   the cut. See `vignette("ddm")` for the measurement.
#' @param allow_unreachable Permit a `max_ndt` above the smallest
#'   response time. Only correct inside a mixture, where another
#'   component carries the rows the Wiener density cannot reach.
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
#' Ratcliff, R. and Tuerlinckx, F. (2002). Estimating parameters of the
#' diffusion model: approaches to dealing with contaminant reaction
#' times and parameter variability. *Psychonomic Bulletin & Review*,
#' 9(3), 438-481.
#'
#' @examples
#' set.seed(1)
#' dat <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
#' fit <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
#'            family = wiener(), data = dat)
#' fixef(fit)
#'
#' @export
wiener <- function(max_ndt = NULL, variability = character(0),
                   nodes = c(sz = 7L, st = 21L),
                   allow_unreachable = FALSE, link = "identity") {
  if (!is.null(max_ndt)) {
    if (!is.numeric(max_ndt) || length(max_ndt) != 1L ||
        !is.finite(max_ndt) || max_ndt <= 0) {
      stop("wiener(): `max_ndt` must be one positive finite number, ",
           "or NULL to take it from the data.", call. = FALSE)
    }
  }
  if (!is.logical(allow_unreachable) || length(allow_unreachable) != 1L ||
      is.na(allow_unreachable)) {
    stop("wiener(): `allow_unreachable` must be TRUE or FALSE.",
         call. = FALSE)
  }
  cfg <- list(max_ndt = max_ndt, link = link,
              allow_unreachable = isTRUE(allow_unreachable),
              variability = ddm_check_variability(variability),
              nodes = ddm_check_nodes(nodes))
  ddm_family(cfg, ub = max_ndt %||% NA_real_, delta = 1e-9)
}

#' The variability names, in the order the dpars are declared.
#'
#' Canonical order rather than the user's, so that two spellings of the
#' same model give the same parameter vector and the same summary.
#'
#' @noRd
ddm_check_variability <- function(variability) {
  known <- c("sv", "sz", "st")
  if (is.null(variability)) variability <- character(0)
  if (!is.character(variability) || anyNA(variability) ||
      !all(variability %in% known) || anyDuplicated(variability)) {
    stop("wiener(): `variability` names the across-trial variability ",
         "parameters to estimate, as a character vector with no ",
         "repeats, drawn from \"sv\" (drift rate), \"sz\" (start ",
         "point) and \"st\" (non-decision time).", call. = FALSE)
  }
  known[known %in% variability]
}

#' Merge a user's node counts over the defaults.
#'
#' @noRd
ddm_check_nodes <- function(nodes) {
  out <- c(sz = 7L, st = 21L)
  if (is.null(nodes) || !length(nodes)) return(out)
  if (!is.numeric(nodes) || is.null(names(nodes)) ||
      !all(names(nodes) %in% names(out)) || anyDuplicated(names(nodes)) ||
      any(!is.finite(nodes)) || any(nodes < 1) ||
      any(nodes != round(nodes))) {
    stop("wiener(): `nodes` gives the Gauss-Legendre node counts as a ",
         "named vector of whole numbers at least 1, with names drawn ",
         "from \"sz\" and \"st\".", call. = FALSE)
  }
  out[names(nodes)] <- as.integer(nodes)
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' A logit scaled onto `(0, ub)`.
#'
#' `ub` may be `NA`, which is the state a family constructed without
#' `max_ndt` is in before [frm()] hands it the data. Using such a link
#' says so rather than returning a silent `NA`.
#'
#' @noRd
ddm_scaled_logit <- function(ub, dpar) {
  force(ub)
  force(dpar)
  bound <- function() {
    if (is.na(ub)) {
      stop("wiener(): the ", dpar, " bound is not set yet. This ",
           "happens when a wiener() family object is used outside ",
           "frm(), for example to inspect its links before a fit. Pass ",
           "`max_ndt` to wiener() to set the bound up front.",
           call. = FALSE)
    }
    ub
  }
  list(
    name = "scaled_logit",
    linkfun = function(mu) { U <- bound(); log(mu / (U - mu)) },
    linkinv = function(eta) bound() / (1 + exp(-eta)),
    mu_eta = function(eta) {
      p <- 1 / (1 + exp(-eta))
      bound() * p * (1 - p)
    })
}

#' Build the family object, given whatever the data has settled.
#'
#' Called twice: once by [wiener()], with the bound the user supplied or
#' `NA`, and once more from `family_finalize()` with the bound and the
#' unreachable-row margin the response determines. The second call is
#' what the fit actually uses.
#'
#' @noRd
ddm_family <- function(cfg, ub, delta) {
  vv <- cfg$variability
  nd <- ddm_nodes(vv, cfg$nodes)
  st_on <- "st" %in% vv
  dpars <- c("mu", "bs", "ndt", "bias", vv)
  links <- list(mu = cfg$link, bs = "log",
                ndt = ddm_scaled_logit(ub, "non-decision time"),
                bias = "logit")
  if ("sv" %in% vv) links$sv <- "log"
  # a width on the same (0, 1) scale as bias, so the logit is the
  # scaled logit its own support asks for
  if ("sz" %in% vv) links$sz <- "logit"
  # and a duration in the response's units, bounded by the same
  # quantity that bounds the non-decision time it is centered on: a log
  # link here lets the optimizer walk out to a width no response time
  # could have come from, where every row's range is cut and the
  # surface is flat
  if (st_on) {
    links$st <- ddm_scaled_logit(if (is.na(ub)) NA_real_ else 2 * ub,
                                 "non-decision time range")
  }

  lpdf <- if (!length(vv) && !cfg$allow_unreachable) {
    # The plain Wiener density, untouched: no variability parameter
    # exists, so there is nothing to average over and no quadrature to
    # pay for, and the bounded link guarantees every row is reachable,
    # so there is nothing to hold off the singularity either.
    function(y, dpars, aterms) {
      ddm_lpdf_both(y - dpars[["ndt"]], dpars[["mu"]], dpars[["bs"]],
                    dpars[["bias"]], ddm_indicator(aterms))
    }
  } else if (!length(vv)) {
    # The same density with the decision time held at `delta`. Only
    # reached when the user has declared that some rows are unreachable,
    # which is the mixture case: there the Wiener component's own
    # likelihood for a fast guess is zero, and it has to be a zero the
    # log-sum-exp can differentiate. -Inf is not: it exponentiates to
    # zero correctly and then contributes NaN to every gradient.
    function(y, dpars, aterms) {
      ddm_lpdf_both(ddm_floor(y - dpars[["ndt"]] - delta, delta),
                    dpars[["mu"]], dpars[["bs"]], dpars[["bias"]],
                    ddm_indicator(aterms))
    }
  } else {
    function(y, dpars, aterms) {
      ddm_lpdf_var(y, dpars[["mu"]], dpars[["bs"]], dpars[["bias"]],
                   dpars[["ndt"]],
                   if (is.null(dpars[["sv"]])) 0 else dpars[["sv"]],
                   if (is.null(dpars[["sz"]])) 0 else dpars[["sz"]],
                   if (is.null(dpars[["st"]])) 0 else dpars[["st"]],
                   ddm_indicator(aterms), nd, st_on, delta)
    }
  }

  init <- list(
    # A drift of zero is the honest starting guess: the sign of the
    # drift is what the data are there to tell us, and a start with
    # the wrong sign costs more than a start at the middle.
    mu = function(y, aterms) 0,
    bs = function(y, aterms) 1.5,
    ndt = function(y, aterms) 0.5 * min(y),
    bias = function(y, aterms) 0.5)
  # Small starts for the variability parameters, because zero is on the
  # boundary of every one of their links and a large start makes the
  # first quadrature straddle a range the data cannot support.
  if ("sv" %in% vv) init$sv <- function(y, aterms) 0.3
  if ("sz" %in% vv) init$sz <- function(y, aterms) 0.05
  if (st_on) init$st <- function(y, aterms) 0.1 * min(y)

  frmtmb::custom_family(
    "wiener",
    dpars = dpars,
    links = links,
    lpdf = lpdf,
    valid_y = function(y, aterms) ddm_check_response(y, aterms),
    family_finalize = function(fam, y, aterms) {
      ddm_finalize(cfg, y)
    },
    init_dpars = init,
    type = "continuous",
    post = list(mean_fn = if (!length(vv)) {
      function(dpars, aterms) ddm_mean_rt(dpars, aterms)
    } else {
      # The closed-form conditional mean is the mean of the WRONG model
      # once the parameters vary between trials, and a post-fit method
      # that quietly returns a number for the wrong model is worse than
      # one that refuses.
      gh <- ddm_gauss_hermite(21L)
      function(dpars, aterms) ddm_mean_rt_var(dpars, aterms, nd, gh)
    }),
    sim = if (!length(vv)) {
      function(dpars, aterms, n) ddm_sim_rt(dpars, aterms, n)
    } else {
      function(dpars, aterms, n) ddm_sim_rt_var(dpars, aterms, n, nd)
    },
    sim_refusal = NULL)
}

#' Fill in everything the family could not know until it had the data.
#'
#' The non-decision-time bound and the margin by which an unreachable
#' row is held off the singularity are both properties of the response.
#' `family_finalize()` is the slot frmtmb provides for exactly this: it
#' runs once at frame assembly, after the response is validated and
#' before any link is used, and whatever it returns is the family the
#' rest of the fit sees. Before that slot existed this package wrote the
#' bound into an environment the link closures read at call time, which
#' worked only for as long as the undocumented call order held.
#'
#' @noRd
ddm_finalize <- function(cfg, y) {
  lo <- min(y)
  ub <- cfg$max_ndt %||% lo
  if (!is.null(cfg$max_ndt) && ub > lo && !cfg$allow_unreachable) {
    # min(y) itself is allowed, and is the default: the scaled logit
    # never reaches its own bound at a finite linear predictor, so
    # ndt < min(rt) stays strict. Anything above min(y) does admit
    # parameter values with no likelihood.
    stop("wiener: max_ndt = ", format(ub), " is above the smallest ",
         "response time (", format(lo), "). The density is zero at ",
         "and below the non-decision time, so a bound above min(rt) ",
         "admits parameter values with no likelihood. In a mixture ",
         "the other component covers those trials, and ",
         "allow_unreachable = TRUE says so.", call. = FALSE)
  }
  ddm_family(cfg, ub = ub, delta = 1e-9 * lo)
}

#' The decision indicator, under whichever spelling supplied it.
#'
#' `dec()` is the term this package contributes to frmtmb's addition-term
#' registry and the one every reference on the model uses; `vint()` is
#' the general-purpose route that was the only one available before that
#' registry existed, and it keeps working.
#'
#' @noRd
ddm_indicator <- function(aterms) {
  up <- aterms[["dec"]]
  if (is.null(up)) up <- aterms[["vint1"]]
  up
}

#' Response and decision-indicator validation.
#'
#' The missing-indicator refusal is still written out here rather than
#' declared through `frmtmb_family(required_aterms =)`, and that is the
#' one hand-rolled check left in this package. `required_aterms` names
#' the terms a density needs ALL of; this family needs EITHER of two,
#' because `dec()` is the spelling to use and `vint()` is the spelling
#' that already works. A declaration cannot say "either".
#'
#' @noRd
ddm_check_response <- function(y, aterms) {
  if (any(!is.finite(y)) || any(y <= 0)) {
    stop("wiener: the response must be a strictly positive, finite ",
         "response time.", call. = FALSE)
  }
  up <- ddm_indicator(aterms)
  if (is.null(up)) {
    # An absent addition term reaches the density as NULL, where the
    # arithmetic silently collapses to a zero-length log-likelihood and
    # the fit "succeeds".
    stop("wiener: the decision indicator is missing. Which boundary a ",
         "trial ended at is data, and it reaches the family through ",
         "dec(), as it does in brms:\n",
         "    frm(bf(rt | dec(decision) ~ x), family = wiener(), ...)\n",
         "where `decision` is a factor whose second level is the upper ",
         "boundary, or a 0/1 column. vint(upper) carries the same ",
         "thing as a plain 0/1 integer and also works.",
         call. = FALSE)
  }
  if (any(!is.finite(up)) || any(up != 0 & up != 1)) {
    stop("wiener: the decision indicator must be 0 (lower boundary) ",
         "or 1 (upper boundary). dec() coerces a factor or a character ",
         "vector for you, taking its second level as the upper ",
         "boundary; vint() does not, so recode it with ",
         "as.integer(decision == \"upper\").", call. = FALSE)
  }
  invisible(NULL)
}
