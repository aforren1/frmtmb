#' The go/no-go diffusion model
#'
#' A two-boundary diffusion in which only ONE boundary is observed. The
#' evidence accumulator runs exactly as in [wiener()], but the
#' experiment records a response only when it reaches the upper
#' boundary, and only if it does so before a deadline. Reaching the
#' lower boundary and never reaching either produce the same
#' observation: no response.
#'
#' This is the design Gomez, Ratcliff and Perea (2007) fitted, and the
#' family EMC2 calls `DDMGNG`. It is the right model whenever one
#' response alternative has no overt response: go/no-go, lexical
#' decision with a single key, signal detection with a single report.
#'
#' Fitting such data with [wiener()] is not an option, because
#' [wiener()] needs to know which boundary each trial ended at and half
#' the trials here do not say.
#'
#' @section The likelihood:
#' A trial contributes one of two things, and which one is data.
#'
#' A GO trial, a response at the upper boundary at time `rt`, contributes
#' the ordinary defective first-passage density of [wiener()] at the
#' upper boundary. It is defective on purpose: it integrates to the
#' probability of a go response, not to one.
#'
#' A NO-GO trial contributes the probability that the upper boundary was
#' not reached before the deadline,
#'
#' \deqn{P(\text{no go}) = 1 - F_{\text{upper}}(TD),}
#'
#' which is the lower-boundary mass accumulated by the deadline plus the
#' mass still diffusing at it. The two add up: integrating the go
#' density from the non-decision time to the deadline gives
#' \eqn{F_{\text{upper}}(TD)}, and adding the no-go probability gives
#' exactly one. Nothing has gone missing, and `test-rdm-gng.R` asserts
#' it.
#'
#' @section Parameters:
#' The same four [wiener()] has, with the same names and the same links,
#' so a model moves between the two families by changing the family and
#' nothing else.
#'
#' \describe{
#'   \item{`mu`}{Drift rate. Identity link, so it is signed: positive
#'     drift favours the go boundary.}
#'   \item{`bs`}{Boundary separation. Log link.}
#'   \item{`ndt`}{Non-decision time. Bounded link, taken from the GO
#'     trials only; see below.}
#'   \item{`bias`}{Relative start point in `(0, 1)`. Logit link. 0.5 is
#'     unbiased.}
#' }
#'
#' The go boundary is the UPPER one, which is a convention rather than a
#' restriction: a model in which the observed response is the lower
#' boundary is this one with the drift negated and the bias reflected.
#'
#' @section The data: which trials, and by when:
#' Two things are per-row data. Which trials produced a response reaches
#' the family through `dec()`, and the deadline reaches it either as one
#' number on the family or through `vreal()`:
#'
#' ```
#' frm(bf(rt | dec(responded) ~ cond), family = wiener_gng(deadline = 1.5),
#'     data = dat)
#'
#' frm(bf(rt | dec(responded) + vreal(deadline) ~ cond),
#'     family = wiener_gng(), data = dat)
#' ```
#'
#' `dec()` here means RESPONDED, coded 1, or did not, coded 0. That is
#' the same 0/1 the term carries for [wiener()] and it means the same
#' thing, because with one observable boundary "reached the upper
#' boundary" and "responded" are the same event. It takes a factor, a
#' character vector or a logical the way brms does, reading the second
#' level as a response, so a plain logical column works as it reads.
#'
#' Give the deadline on the family when every trial had the same one,
#' which is the usual design, and through `vreal()` when it varied.
#' Supplying neither is refused by name.
#'
#' A no-go trial has no response time. The response column still needs a
#' number, because it is the response and frmtmb will not carry an `NA`
#' through a rowwise family; the likelihood does not read it, and the
#' deadline is the natural thing to put there. Any positive finite value
#' gives the same answer.
#'
#' A go trial whose response time is past its own deadline is refused
#' rather than fitted, because the model gives that trial no
#' probability at all: the deadline is what stopped the accumulator.
#'
#' @section Non-decision time:
#' The bound is the fastest GO response, not the fastest row. A no-go
#' row's response entry is a placeholder, so letting it into the bound
#' would let a placeholder decide a parameter's range. A model with no
#' go trials at all is refused: nothing in it identifies the
#' non-decision time.
#'
#' @section Across-trial variability, and why there is none here:
#' [wiener()] offers Ratcliff's `sv`, `sz` and `st0`. This family offers
#' none of the three, and the reason is the no-go branch rather than an
#' oversight.
#'
#' The go branch would inherit all three for free, because it is
#' [wiener()]'s density. The no-go branch would not. The drift enters
#' the DENSITY as an exponential-quadratic, which is what makes
#' averaging over a normal drift exact and free; it enters the
#' DISTRIBUTION FUNCTION through the eigenvalues of both series, where
#' it is not, so `sv` would need a quadrature of its own rather than a
#' completed square. `sz` and `st0` could be reached with the existing
#' Gauss-Legendre nodes, but shipping two of three would make
#' `variability =` mean something different here than it does on
#' [wiener()], which is worse than not having it.
#'
#' EMC2's `DDMGNG` does carry all three, by calling a compiled
#' distribution function that integrates them numerically. If you need
#' them, that is where they are.
#'
#' @section Accuracy:
#' The no-go branch needs the Wiener defective distribution function,
#' which this package did not have: [wiener()] declares no `lcdf`, and
#' says so in its own compatibility table. It is written for this family
#' in `R/wiener-cdf.R`, as two series blended in `log(u)` the way the
#' density's two are.
#'
#' The blend is centred much lower than the density's, at `u = 0.02`
#' rather than 0.35, and that is the substantive choice here. The
#' small-time route reaches the no-go probability as `1 - F_upper`,
#' which cancels precisely where a no-go trial is surprising; the
#' large-time route computes it directly, as a gambler's-ruin
#' probability plus a correction, and never subtracts. Handing over as
#' early as the large-time route is accurate is what keeps the
#' likelihood honest on the rows that pull hardest on it.
#'
#' Measured against a 260-bit reference over 1200 points spanning
#' `t` in 0.05 to 15, drift in -2 to 5, boundary separation 0.8 to 4
#' and relative start point 0.25 to 0.9, the blended probability holds
#' **2.2e-12** relative, with one row worse than 1e-12 and none worse
#' than 1e-9.
#'
#' The centre was 0.06 at 0.3.0 and is 0.02 because of what a wider
#' grid showed. At `t = 2.5`, drift 5, separation 4 and start point
#' 0.90, where the true no-go probability is 8.2e-16, the 0.06 blend
#' was 7.75e-05 relative. Neither series was at fault: the large-time
#' route alone was 6.7e-15 there, and the small-time route was 12.1
#' RELATIVE, because `1 - F_upper` had nothing left to subtract from.
#' A weight of about 3e-05 on a hopeless number is what the error was.
#' A smooth blend has no safe side unless both branches degrade
#' gently, and this one has a branch that does not, so the weight has
#' to saturate before that branch collapses.
#'
#' The two series are independent derivations and agree with each other
#' to 4.8e-78 at 260 bits, which is the check that the derivation is
#' right rather than merely stable. Against `WienR::pWDM()`, which is
#' what EMC2 calls, agreement is 3.3e-13 wherever the no-go probability
#' is above 0.01. Below that the reference is the weaker of the two,
#' though by how much depends on what it is asked for. At one grid
#' point with a no-go probability of 1.19e-13, `1 - WienR::pWDM()` is
#' 4.4 percent wrong at WienR's DEFAULT precision and at every setting
#' down to `precision = 1e-12`, 0.13 percent at 1e-14 and 0.033 percent
#' at 1e-16. What matters for the comparison with EMC2 is the first
#' row of that list: `EMC2:::pDDM` calls `pWDM` with
#' `precision = 0.005`, looser than any of them, so EMC2 in practice
#' sits at the 4.4 percent end. This family's large-time route is
#' 3.0e-15 there.
#'
#' @param deadline The response deadline, in the units of the response.
#'   One positive number when every trial shares it. `NULL`, the
#'   default, takes it per row from `vreal()`.
#' @param max_ndt Upper bound for the non-decision time. `NULL`, the
#'   default, takes it from the fastest go response.
#'
#' @return A `frmtmb_family`.
#'
#' @references
#' Gomez, P., Ratcliff, R. and Perea, M. (2007). A model of the go/no-go
#' task. *Journal of Experimental Psychology: General*, 136(3), 389-413.
#'
#' @seealso [wiener_gng_simulate()] to generate from the model, and
#'   [wiener()] for the two-choice diffusion this is the one-boundary
#'   version of.
#'
#' @examples
#' set.seed(1)
#' dat <- wiener_gng_simulate(600, mu = 1.0, bs = 1.4, ndt = 0.25,
#'                            deadline = 1.5)
#' mean(dat$responded)
#' fit <- frm(bf(rt | dec(responded) ~ 1, bias = 0.5),
#'            family = wiener_gng(deadline = 1.5), data = dat)
#' fixef(fit)
#'
#' @export
wiener_gng <- function(deadline = NULL, max_ndt = NULL) {
  if (!is.null(deadline)) {
    if (!is.numeric(deadline) || length(deadline) != 1L ||
        !is.finite(deadline) || deadline <= 0) {
      stop("wiener_gng(): `deadline` is the time at which a trial stops ",
           "waiting for a response and must be one positive finite ",
           "number, or NULL to take it per row from vreal().",
           call. = FALSE)
    }
    deadline <- as.numeric(deadline)
  }
  if (!is.null(max_ndt)) {
    if (!is.numeric(max_ndt) || length(max_ndt) != 1L ||
        !is.finite(max_ndt) || max_ndt <= 0) {
      stop("wiener_gng(): `max_ndt` bounds the non-decision time and ",
           "must be one positive finite number, or NULL to read it off ",
           "the go responses.", call. = FALSE)
    }
  }

  need <- if (is.null(deadline)) c("dec", "vreal1") else "dec"

  fam <- frmtmb::custom_family(
    "wiener_gng",
    dpars = c("mu", "bs", "ndt", "bias"),
    links = list(mu = "identity", bs = "log", ndt = "log",
                 bias = "logit"),
    lpdf = function(y, dpars, aterms) {
      gng_lpdf(y, dpars, aterms, deadline)
    },
    valid_y = function(y, aterms) {
      gng_check_response(y, aterms, deadline)
    },
    family_finalize = function(fam, y, aterms) {
      gng_finalize(fam, y, aterms, max_ndt)
    },
    required_aterms = need,
    init_dpars = list(
      mu = function(y, aterms) 1,
      bs = function(y, aterms) log(1.4),
      # The bound is not known here, so the start is expressed against
      # the go responses the bound will be taken from.
      ndt = function(y, aterms) 0.5 * min(gng_go_rt(y, aterms)),
      bias = function(y, aterms) 0.5),
    type = "continuous",
    # A refusing mean rather than no mean at all, for the reason rdm()
    # gives and with a sharper version of it here: this family's primary
    # dpar IS called `mu`, so with the slot empty frmtmb returned the
    # DRIFT RATE from fitted() and predict(type = "response"), a
    # constant 1.05 for data whose go response times average 0.6. The
    # refusal is also the honest answer, not merely a safer one; see
    # below.
    post = list(mean_fn = function(dpars, aterms) {
      stop("wiener_gng: this family has no mean response to report. A ",
           "go/no-go trial produces a PAIR, whether a response ",
           "happened and when, and no single number summarizes it: the ",
           "no-go trials have no response time to average at all. ",
           "fitted(), predict(type = \"response\") and residuals(type ",
           "= \"response\") are unavailable for that reason rather ",
           "than for want of an integral. predict(type = \"link\") ",
           "gives the drift, boundary separation, non-decision time ",
           "and bias the fit estimated.", call. = FALSE)
    }),
    sim = function(dpars, aterms, n) {
      gng_sim_rt(dpars, aterms, n, deadline)
    },
    primary_dpars = "mu")

  # Carried so that a reader of the fitted object can see the deadline
  # the likelihood used, and so the simulator need not be told again.
  fam[["gng_deadline"]] <- deadline
  fam
}

#' The two branches, selected by a 0/1 mask rather than a branch.
#'
#' The winner of a race and the responder of a go/no-go trial are both
#' data, so both families select the same way: multiply by an indicator
#' and add. One tape serves every row.
#'
#' Both times are floored before use. A no-go row's response entry is a
#' placeholder the go branch still evaluates at, and a `-Inf` there
#' would meet the mask's zero and give `NaN`, which is the one thing a
#' zero weight cannot absorb.
#'
#' @noRd
gng_lpdf <- function(y, dpars, aterms, deadline) {
  go <- aterms[["dec"]]
  td <- if (is.null(deadline)) aterms[["vreal1"]] else deadline
  ndt <- dpars[["ndt"]]
  lgo <- ddm_lpdf_both(ddm_floor(y - ndt, 1e-12), dpars[["mu"]],
                       dpars[["bs"]], dpars[["bias"]], 1)
  lng <- ddm_nogo_lprob(ddm_floor(td - ndt, 1e-12), dpars[["mu"]],
                        dpars[["bs"]], dpars[["bias"]])
  go * lgo + (1 - go) * lng
}

#' The response times of the go trials alone.
#'
#' The bound on the non-decision time and its starting value are both
#' properties of the trials that produced a response. A no-go row's
#' entry is a placeholder, and letting a placeholder into `min()` would
#' let it decide a parameter's range.
#'
#' @noRd
gng_go_rt <- function(y, aterms) {
  go <- as.numeric(aterms[["dec"]])
  y[!is.na(go) & go > 0.5]
}

#' Response, indicator and deadline validation.
#'
#' @noRd
gng_check_response <- function(y, aterms, deadline) {
  if (any(!is.finite(y)) || any(y <= 0)) {
    stop("wiener_gng: the response must be strictly positive and ",
         "finite on every row. A no-go trial has no response time, but ",
         "the column still needs a number the likelihood will ignore; ",
         "the deadline is the natural one to put there.", call. = FALSE)
  }
  if (!is.null(aterms[["vint1"]])) {
    stop("wiener_gng: the response indicator travels through dec(), ",
         "not vint(). dec() says whether a trial produced a response, ",
         "coded 1, or did not, coded 0, and it reads a factor or a ",
         "logical column the way brms does. vint() is free for ",
         "anything else the model needs.", call. = FALSE)
  }
  go <- aterms[["dec"]]
  if (is.null(go) || any(!is.finite(go)) || any(!go %in% c(0, 1))) {
    stop("wiener_gng: dec() must be 0 or 1 on every row, 1 for a trial ",
         "that produced a response. A factor, a character vector or a ",
         "logical is coerced on its levels with the second one read as ",
         "a response.", call. = FALSE)
  }
  if (!any(go == 1)) {
    stop("wiener_gng: no trial in these data produced a response, so ",
         "nothing identifies the non-decision time or the shape of the ",
         "response time distribution. A go/no-go model needs go ",
         "trials.", call. = FALSE)
  }
  td <- if (is.null(deadline)) aterms[["vreal1"]] else deadline
  if (any(!is.finite(td)) || any(td <= 0)) {
    stop("wiener_gng: the deadline must be positive and finite on ",
         "every row. It is the time at which a trial stops waiting, so ",
         "a deadline of zero or less leaves no window in which a ",
         "response could have happened.", call. = FALSE)
  }
  late <- go == 1 & y > rep_len(as.numeric(td), length(y)) *
    (1 + 1e-9)
  if (any(late)) {
    stop("wiener_gng: ", sum(late), " trial(s) recorded a response ",
         "after their own deadline, the fastest at ",
         format(min(y[late])), " against a deadline of ",
         format(min(rep_len(as.numeric(td), length(y))[late])),
         ". The deadline is what stops the accumulator, so the model ",
         "gives those trials no probability at all. Check that the ",
         "deadline is in the units of the response.", call. = FALSE)
  }
  invisible(NULL)
}

#' Fit the non-decision-time link to the GO responses.
#'
#' @noRd
gng_finalize <- function(fam, y, aterms, max_ndt) {
  ddm_ndt_finalize(fam, gng_go_rt(y, aterms), max_ndt, "wiener_gng")
}

#' Draw a response time conditional on each row's own outcome.
#'
#' A no-go row has no response time to draw, so it gets its deadline,
#' which is the placeholder the family documents. A go row is drawn by
#' rejection from the unconditional process until it produces a
#' response before its deadline, which is exactly the go density
#' renormalized on its own support.
#'
#' @noRd
gng_sim_rt <- function(dpars, aterms, n, deadline) {
  go <- rep_len(as.numeric(aterms[["dec"]]), n)
  td <- rep_len(as.numeric(
    if (is.null(deadline)) aterms[["vreal1"]] else deadline), n)
  mu <- rep_len(as.numeric(dpars[["mu"]]), n)
  bs <- rep_len(as.numeric(dpars[["bs"]]), n)
  ndt <- rep_len(as.numeric(dpars[["ndt"]]), n)
  bias <- rep_len(as.numeric(dpars[["bias"]]), n)
  out <- td
  todo <- which(go > 0.5)
  for (round in 1:200) {
    if (!length(todo)) break
    d <- ddm_simulate(length(todo), mu = mu[todo], bs = bs[todo],
                      ndt = ndt[todo], bias = bias[todo])
    hit <- d$upper == 1 & d$rt <= td[todo]
    out[todo[hit]] <- d$rt[hit]
    todo <- todo[!hit]
  }
  if (length(todo)) {
    stop("wiener_gng: could not draw a response time for ",
         length(todo), " go trial(s). The fitted parameters make a ",
         "response before the deadline so unlikely that 200 rounds of ",
         "rejection produced none. Use wiener_gng_simulate() for an ",
         "unconditional draw, which does not have to hit a given ",
         "outcome.", call. = FALSE)
  }
  out
}

#' Simulate from a go/no-go diffusion model
#'
#' Draws from the generative process: run the two-boundary diffusion,
#' and record a response only if it reached the upper boundary before
#' the deadline. Nothing here evaluates the go density or the no-go
#' probability, so the simulator is an independent statement of the same
#' model and can be used to check the likelihood rather than merely to
#' agree with it.
#'
#' @param n Number of trials.
#' @param mu Drift rate.
#' @param bs Boundary separation.
#' @param ndt Non-decision time.
#' @param bias Relative start point in `(0, 1)`.
#' @param deadline The response deadline. One number, or one per trial.
#'
#' @return A data frame with `responded` (1 for a go trial), `rt` (the
#'   response time on a go trial, and the deadline on a no-go trial,
#'   which is the placeholder [wiener_gng()] documents) and `deadline`.
#'
#' @examples
#' set.seed(1)
#' dat <- wiener_gng_simulate(1000, mu = 1.0, bs = 1.4, ndt = 0.25,
#'                            deadline = 1.5)
#' mean(dat$responded)
#' summary(dat$rt[dat$responded == 1])
#'
#' @export
wiener_gng_simulate <- function(n, mu = 1, bs = 1.4, ndt = 0.25,
                                bias = 0.5, deadline = 1.5) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("wiener_gng_simulate(): `n` is the number of trials to draw, ",
         "one whole number of 1 or more.", call. = FALSE)
  }
  n <- as.integer(n)
  td <- rep_len(as.numeric(deadline), n)
  if (any(!is.finite(td)) || any(td <= 0)) {
    stop("wiener_gng_simulate(): `deadline` must be positive and ",
         "finite, as one number or one per trial.", call. = FALSE)
  }
  d <- ddm_simulate(n, mu = mu, bs = bs, ndt = ndt, bias = bias)
  responded <- d$upper == 1 & d$rt <= td
  data.frame(responded = as.integer(responded),
             rt = ifelse(responded, d$rt, td),
             deadline = td)
}
