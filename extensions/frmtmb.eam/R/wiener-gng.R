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
#' @section Across-trial variability:
#' Ratcliff's `sv`, `sz` and `st` are here, named and linked exactly as
#' [wiener()] names and links them, so a model moves between the two
#' families by changing the family and nothing else:
#'
#' ```
#' frm(bf(rt | dec(responded) ~ cond, bias = 0.5),
#'     family = wiener_gng(deadline = 1.5,
#'                         variability = c("sv", "sz", "st")),
#'     data = dat)
#' ```
#'
#' The GO branch is [wiener()]'s density and inherits all three
#' unchanged: it calls the same averaged density at the same nodes, so
#' the two families agree to the last bit on a go trial. The NO-GO
#' branch is the part that had to be written, and it is a different
#' integral, because it averages the DISTRIBUTION FUNCTION rather than
#' the density.
#'
#' Each of the three enters it its own way.
#'
#' \describe{
#'   \item{`sv`}{Quadrature, and the only one of the three that costs
#'     more here than in [wiener()]. The drift enters the density only
#'     as `exp(-v a w - v^2 t / 2)`, an exponential-quadratic that a
#'     normal average integrates in closed form. It enters the
#'     distribution function through the eigenvalues as well, as
#'     `1 / (v^2 a^2 + k^2 pi^2)` in every term of the large-time
#'     series and as the gambler's-ruin probability that series
#'     corrects, and neither is an exponential-quadratic. So the no-go
#'     branch integrates the drift by Gauss-Hermite, and `nogo_nodes`
#'     carries an `sv` entry that [wiener()] has no use for.}
#'   \item{`sz`}{The same Gauss-Legendre nodes the density uses, over
#'     the same uniform start point.}
#'   \item{`st`}{Shifts the DEADLINE. The density's non-decision-time
#'     range is cut at the response time, because past that the decision
#'     time is negative and the density is zero; the no-go probability
#'     has no such cut, because a non-decision time past the deadline
#'     leaves the accumulator no time at all and the probability of no
#'     response is then exactly one. The whole range is averaged.}
#' }
#'
#' The identity that ties the two branches together holds WHILE THE
#' START-POINT RANGE STAYS INSIDE THE BOUNDARIES, which is where the
#' model is defined: the go density integrated to the deadline plus the
#' no-go probability is one, to 8.9e-16 in the plain family and to
#' whatever the quadrature gives once a variability parameter is on.
#' With the default `nogo_nodes` that is 3.9e-14 at `sv` = 0.6 and
#' 2.0e-10 at `st` = 0.20, but only 8.1e-05 at `sv` = 2.0 and 9.9e-08
#' at `st` = 0.45. It is not "exactly" one, and the size of the gap is
#' the `sv` row of the node table below.
#'
#' Outside the boundaries the two branches are no longer the same
#' average of the same pair, and the mass is not conserved. Only the
#' NO-GO branch clamps the start point; the go branch is left unclamped
#' deliberately, because clamping it would break the bit-identity with
#' [wiener()]. Measured, at a deadline of 1.5:
#'
#' | | go | no-go | total |
#' |---|---|---|---|
#' | `sz` = 0.5, `bias` = 0.85 | 0.919017 | 0.051622 | 0.970639 |
#' | `sz` = 0.9, `bias` = 0.90 | 0.780184 | 0.063186 | 0.843370 |
#'
#' so up to 16 percent of the pair's mass is missing there. What saves
#' it is that the region is strictly downhill. Profiling `sz` on 2000
#' trials generated at `bias` = 0.85 and a true `sz` of 0.20, where the
#' range leaves the boundary at `sz` = 0.30, the best point inside is
#' 1867.975 at `sz` = 0.190 and the best point outside is 1678.640 at
#' `sz` = 0.305: the interior peak wins by 189 log units and the surface
#' is monotone across the crossing. The clamp is a barrier an optimizer
#' walks away from, not an attractor, so a fit does not end up there.
#' Read a fitted `sz` whose range crosses a boundary as a fit that has
#' run out of model, not as an estimate.
#'
#' @section How many nodes, and what they cost:
#' Measured against a 200-bit reference, over six parameter settings
#' spanning deadlines 0.8 to 2.5, drifts 0.5 to 3, boundary separations
#' 0.8 to 2.5 and start points 0.4 to 0.8. Worst relative error, one
#' parameter live at a time:
#'
#' | nodes | `sz` = 0.1 | `sz` = 0.3 | | `st` = 0.05 | `st` = 0.3 |
#' |---|---|---|---|---|---|
#' | 3 | 4.5e-06 | 2.1e-03 | | 3.4e-11 | 1.6e-06 |
#' | 5 | 1.8e-11 | 6.4e-07 | | 8.6e-15 | 4.0e-10 |
#' | 7 | 5.1e-15 | 4.2e-11 | | 6.3e-15 | 6.1e-14 |
#' | 9 | 4.0e-15 | 2.1e-15 | | 4.2e-15 | 2.7e-15 |
#'
#' Both saturate by seven nodes, and `st` saturates there where the
#' DENSITY's own `st` integral needs 21. That is not a discrepancy: the
#' density's range is cut at the response time and its integrand turns
#' on sharply at the cut, and the probability has no cut at all.
#'
#' `sv` is the one that does not saturate:
#'
#' | nodes | `sv` = 0.3 | `sv` = 0.8 | `sv` = 1.5 |
#' |---|---|---|---|
#' | 7 | 9.0e-07 | 8.9e-04 | 1.8e-02 |
#' | 11 | 9.5e-12 | 3.8e-06 | 1.0e-02 |
#' | 15 | 8.4e-15 | 1.3e-06 | 1.4e-03 |
#' | 21 | 9.5e-15 | 3.3e-09 | 3.5e-05 |
#' | 31 | 1.2e-14 | 1.2e-12 | 6.6e-07 |
#' | 41 | 1.4e-14 | 2.8e-14 | 9.1e-08 |
#'
#' The node count a given accuracy needs rises with the product of the
#' boundary separation and `sv`, because that product sets how sharp the
#' transition in the drift is. The default of 15 holds 1e-14 at a narrow
#' drift variability and 1e-6 at a moderate one; **raise
#' `nogo_nodes = c(sv = 31)` or higher for a wide one**, and read the
#' table rather than assuming the default is enough.
#'
#' The three do not compound. Measured jointly, the error of the full
#' three-dimensional rule tracks the `sv` error alone to within a factor
#' of three at every setting, so `sz` and `st` at seven nodes are there
#' to not be the binding term, and `sv` is the only knob worth turning.
#'
#' The cost is a product, and it is the reason `nogo_nodes` exists as an
#' argument separate from `nodes`. On 500 rows with all three live:
#'
#' | `nogo_nodes` | grid | wall |
#' |---|---|---|
#' | 7, 5, 5 | 175 | 124 s |
#' | 11, 5, 5 | 275 | 185 s |
#' | 11, 7, 7 | 539 | 610 s |
#'
#' A model with one variability parameter pays one dimension of that and
#' is cheap; the three-at-once model is the expensive one, and it is
#' expensive because a distribution function with no closed form in any
#' of the three has to be evaluated on a product grid.
#'
#' EMC2's `DDMGNG` carries the same three parameters, by calling a
#' compiled distribution function that integrates them numerically. The
#' two agree: over 108 grid points the no-go probability matches
#' `1 - EMC2:::pDDM()` to 9.4e-14 with no variability, 2.4e-13 under
#' `sv` and 3.9e-13 under `sz`, and the go density matches
#' `EMC2:::dDDM()` to 1.2e-15.
#'
#' `st` needs a shift before the two are comparable, and the difference
#' is a CONVENTION rather than a defect on either side. EMC2 takes
#' `st0` as a uniform on `[t0, t0 + st0]`; this family takes `st`
#' centred on `ndt`, following brms. Compare EMC2 at `t0` with this
#' family at `ndt = t0 + st0 / 2` and the agreement is 1.5e-13 for the
#' probability and 1.6e-15 for the density. Compare them without the
#' shift and they differ by 14 to 45 percent, in the density as much as
#' in the probability.
#'
#' @section What a go/no-go design can and cannot identify:
#' `sv` is weakly identified here, and it is the design rather than the
#' likelihood. With only one boundary observed, drift variability trades
#' against the drift and the boundary separation along a ridge.
#'
#' Measured on 3000 simulated trials with a true `sv` of 0.6, this
#' family returns 0.001 while [wiener()] on the same generative
#' parameters, seeing BOTH boundaries, returns 0.612. Profiled at the
#' true values of everything else the likelihood does peak at the truth,
#' but the fit reaches a HIGHER log likelihood at `sv` near zero by
#' moving the drift from 1.20 to 1.09 and the separation from 1.40 to
#' 1.33. Nothing is wrong with the surface; the sample simply prefers
#' that corner of the ridge.
#'
#' How much of that is the ONE observable boundary is less settled than
#' the contrast suggests, and the honest statement is narrower. A grid
#' search on two-boundary data from the same generative parameters, with
#' `ndt` and `bias` held at the truth, finds the same corner there too:
#' `sv` near zero beats the truth by 2.4 log units, against 1.8 on the
#' go/no-go data. So the ridge is a property of the drift-diffusion
#' likelihood and not only of this design; what the go/no-go design does
#' is remove the boundary-proportion information that would otherwise
#' help pin `sv`, and the evidence for that is the pair of fits above
#' rather than a study. Either way the practical advice is the same.
#'
#' So read a small fitted `sv` here as "the data did not pin it" rather
#' than as "there is no drift variability", and prefer to fix it, or to
#' estimate it from a two-choice condition of the same experiment, over
#' reading it off a go/no-go block. `sz` and `st` are better behaved,
#' because both change the SHAPE of the go response time distribution
#' rather than trading against its location.
#'
#' @section Censoring:
#' A trial whose clock was stopped before it responded is a
#' RIGHT-CENSORED observation, and the probability of it is the no-go
#' probability at the time the clock stopped. This family declares that
#' as its log survivor function, so `cens()` works:
#'
#' ```
#' frm(bf(rt | dec(responded) + cens(stopped) ~ cond),
#'     family = wiener_gng(deadline = 1.5), data = dat)
#' ```
#'
#' Left censoring, interval censoring and `trunc()` are refused, and not
#' for want of an integral. The likelihood here is a defective density
#' plus a point mass at "no response", and a truncation window on the
#' response scale renormalizes the density while saying nothing about
#' the mass, so the two halves of every row would be divided by
#' different things. Right censoring has no such problem, because it
#' replaces a whole row rather than reweighting it. So this family
#' declares an `lccdf` and, deliberately, no `lcdf`.
#'
#' The refusal you will see is frmtmb's own and does not name this
#' family:
#'
#' ```
#' cens()/trunc() need a family with a CDF (currently: gaussian,
#' lognormal, poisson, exponential, weibull, inverse.gaussian, cox) ...
#' ```
#'
#' Read it as a decision rather than as an omission. It arrives from
#' frame assembly, which runs before any family-supplied check, so this
#' family has no seam that does not require claiming a CDF it does not
#' have.
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
#' @param variability Which across-trial variability parameters to
#'   estimate: any of `"sv"` (drift rate), `"sz"` (start point) and
#'   `"st"` (non-decision time), exactly as [wiener()] takes them. The
#'   default estimates none.
#' @param nodes Gauss-Legendre node counts for the GO branch, which is
#'   [wiener()]'s density: the same names, the same meaning and the same
#'   defaults, so a go trial is scored identically by the two families.
#' @param nogo_nodes Quadrature node counts for the NO-GO probability,
#'   which is a different integral over the same three distributions.
#'   `sz` and `st` are Gauss-Legendre counts, and `st` defaults lower
#'   than the density's because the probability's range is not cut by
#'   the response time; `sv` is a Gauss-Hermite count with no
#'   counterpart in [wiener()], because there the drift integral is
#'   closed form. Only the entries for the `variability` parameters in
#'   use are read. Raise `sv` for a wide drift variability; see
#'   Across-trial variability for the measured table.
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
wiener_gng <- function(deadline = NULL, max_ndt = NULL,
                       variability = character(0),
                       nodes = c(sz = 7L, st = 21L),
                       nogo_nodes = c(sv = 15L, sz = 7L, st = 7L)) {
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
  cfg <- list(deadline = deadline, max_ndt = max_ndt,
              variability = ddm_check_variability(variability,
                                                  "wiener_gng"),
              nodes = ddm_check_nodes(nodes),
              nogo_nodes = gng_check_nodes(nogo_nodes))
  gng_family(cfg, ub = max_ndt %||% NA_real_, delta = 1e-9)
}

#' Merge a user's node counts over the NO-GO branch's defaults.
#'
#' Why the no-go branch has its own counts at all, when it averages the
#' same three distributions the density averages: they are two different
#' integrands, and forcing them to share a rule is what makes the model
#' unfittable.
#'
#' The density's `st` count is 21 because ITS range is cut at the
#' response time and the integrand turns on sharply at the cut. The
#' probability has no cut, so seven nodes reach 6.1e-14 where the
#' density needs 21 to reach anything comparable. That would be a
#' curiosity if the two branches cost the same, and it is not: the
#' density's rule is a two-dimensional product and the probability's is
#' three-dimensional, so an `st` count carried over from the density is
#' multiplied by every other node count in the grid. Measured, at 500
#' rows, the density's own counts here give 1617 evaluations of a
#' two-series distribution function per row and the tape runs out of
#' memory.
#'
#' `sv` has no counterpart in [wiener()] at all. The drift enters the
#' density only as an exponential-quadratic, which a normal average
#' integrates in closed form; it enters the distribution function
#' through the eigenvalues as well, where nothing does. It is the one
#' count here that buys accuracy per node rather than saturating, and
#' the accuracy it buys falls off as the product of the boundary
#' separation and `sv` grows. See `?wiener_gng` for the measured table
#' and for what to raise it to.
#'
#' @noRd
gng_check_nodes <- function(nodes) {
  out <- c(sv = 15L, sz = 7L, st = 7L)
  if (is.null(nodes) || !length(nodes)) return(out)
  if (!is.numeric(nodes) || is.null(names(nodes)) ||
      !all(names(nodes) %in% names(out)) || anyDuplicated(names(nodes)) ||
      any(!is.finite(nodes)) || any(nodes < 1) ||
      any(nodes != round(nodes))) {
    stop("wiener_gng(): `nogo_nodes` gives the no-go probability's ",
         "quadrature node counts as a named vector of whole numbers at ",
         "least 1, with names drawn from \"sv\", \"sz\" and \"st\".",
         call. = FALSE)
  }
  out[names(nodes)] <- as.integer(nodes)
  out
}

#' The node set the NO-GO branch integrates on.
#'
#' The go branch gets `ddm_nodes(variability, cfg$nodes)`, which is
#' exactly what [wiener()] gets, so the go density is the same function
#' at the same points. This one is the other rule.
#'
#' An integral that is switched off gets the ONE-node rule rather than
#' no rule, and the family passes it a width of zero, so the switched-off
#' case is not an approximation of anything: one node at the middle of a
#' zero-width interval is an evaluation at the middle, exactly. That is
#' what makes the no-variability family bit-identical to 0.3.0 through
#' the same code path.
#'
#' @noRd
gng_nodes <- function(variability, nodes) {
  nd <- ddm_nodes(variability, nodes)
  nd$sv <- if ("sv" %in% variability) {
    ddm_gauss_hermite(nodes[["sv"]])
  } else {
    list(x = 0, w = 1)
  }
  nd
}

#' Build the family object, given whatever the data has settled.
#'
#' Called twice, as [wiener()]'s builder is: once by [wiener_gng()] with
#' whatever the user supplied, and once from `family_finalize()` with
#' the non-decision-time bound and the unreachable-row margin the go
#' responses determine. The second call is what the fit uses.
#'
#' `delta` matters only under `variability`, and it is [wiener()]'s
#' `delta` for a reason: it is the margin by which the non-decision-time
#' range is held below each row's own response time, so a go trial
#' scored here and the same trial scored by [wiener()] agree to the last
#' bit only if the two margins agree.
#'
#' @noRd
gng_family <- function(cfg, ub, delta) {
  deadline <- cfg$deadline
  vv <- cfg$variability
  # Two rules, because the two branches are two integrands. The go
  # branch's is wiener()'s, name for name and count for count.
  nd <- ddm_nodes(vv, cfg$nodes)
  ndc <- gng_nodes(vv, cfg$nogo_nodes)
  st_on <- "st" %in% vv
  dpars <- c("mu", "bs", "ndt", "bias", vv)
  links <- list(mu = "identity", bs = "log", ndt = "log", bias = "logit")
  if ("sv" %in% vv) links$sv <- "log"
  if ("sz" %in% vv) links$sz <- "logit"
  if (st_on) {
    links$st <- ddm_scaled_logit(if (is.na(ub)) NA_real_ else 2 * ub,
                                 "non-decision time range")
  }

  need <- if (is.null(deadline)) c("dec", "vreal1") else "dec"

  lpdf <- if (!length(vv)) {
    function(y, dpars, aterms) gng_lpdf(y, dpars, aterms, deadline)
  } else {
    function(y, dpars, aterms) {
      gng_lpdf_var(y, dpars, aterms, deadline, nd, ndc, st_on, delta)
    }
  }

  init <- list(
    mu = function(y, aterms) 1,
    # RESPONSE scale, which is what init_dpars takes: frmtmb puts each
    # value through the dpar's own linkfun (?frmtmb_family, step 4).
    # This read log(1.4) until 0.3.1 and so started the separation at
    # 0.336, a factor of four low. Measured over seven designs including
    # bs = 5.0 and variability-on models, every fit still reached the
    # same optimum; the cost was 0 to 5 extra iterations. A wrong start
    # against a documented contract is worth one word whatever it costs.
    bs = function(y, aterms) 1.4,
    # The bound is not known here, so the start is expressed against
    # the go responses the bound will be taken from.
    ndt = function(y, aterms) 0.5 * min(gng_go_rt(y, aterms)),
    bias = function(y, aterms) 0.5)
  # Small starts, for the reason wiener() gives: zero is on the boundary
  # of every one of these links, and a large start makes the first
  # quadrature straddle a range the data cannot support.
  if ("sv" %in% vv) init$sv <- function(y, aterms) 0.3
  if ("sz" %in% vv) init$sz <- function(y, aterms) 0.05
  if (st_on) init$st <- function(y, aterms) 0.1 * min(gng_go_rt(y, aterms))

  fam <- frmtmb::custom_family(
    "wiener_gng",
    dpars = dpars,
    links = links,
    lpdf = lpdf,
    valid_y = function(y, aterms) {
      gng_check_response(y, aterms, deadline)
    },
    family_finalize = function(fam, y, aterms) {
      gng_finalize(cfg, y, aterms)
    },
    required_aterms = need,
    init_dpars = init,
    type = "continuous",
    # Right censoring only, and by declaration. The log survivor is this
    # family's own no-go probability with the deadline moved to the
    # censoring time, which is exactly what "had not responded by then"
    # means, so the seam is a rename rather than a derivation. `lcdf` is
    # deliberately absent, and frmtmb then refuses left censoring,
    # interval censoring and trunc() BY NAME. They are not merely
    # unwritten: this likelihood is a defective density plus a point
    # mass, and a window normalizer computed on the response scale alone
    # would renormalize the no-go rows against a window they do not live
    # in.
    lccdf = function(q, dpars, aterms) gng_lccdf(q, dpars, ndc),
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

#' The same two branches with the across-trial variability integrated
#' out.
#'
#' The go branch is [wiener()]'s averaged density, called with `up = 1`
#' and with the same node sets and the same `delta`, so a go trial
#' scored here and the same trial scored by [wiener()] are the same
#' arithmetic in the same order. The no-go branch is the new integral.
#'
#' A no-go row's response entry is still a placeholder the go branch
#' evaluates at, and it is still ignored: `ddm_lpdf_var()` holds every
#' node `delta` below that row's own response time, so the go branch
#' returns a finite number there whatever the placeholder is, and the
#' mask multiplies it by zero.
#'
#' @noRd
gng_lpdf_var <- function(y, dpars, aterms, deadline, nd, ndc, st_on,
                         delta) {
  go <- aterms[["dec"]]
  td <- if (is.null(deadline)) aterms[["vreal1"]] else deadline
  gv <- function(nm) if (is.null(dpars[[nm]])) 0 else dpars[[nm]]
  sv <- gv("sv")
  sz <- gv("sz")
  st <- gv("st")
  lgo <- ddm_lpdf_var(y, dpars[["mu"]], dpars[["bs"]], dpars[["bias"]],
                      dpars[["ndt"]], sv, sz, st, 1, nd, st_on, delta)
  lng <- ddm_nogo_lprob_var(td, dpars[["mu"]], dpars[["bs"]],
                            dpars[["bias"]], dpars[["ndt"]],
                            sv, sz, st, ndc)
  go * lgo + (1 - go) * lng
}

#' Log probability that no response had happened by time `q`.
#'
#' The no-go probability, with the deadline replaced by the censoring
#' time. Nothing is derived for this: a right-censored go/no-go trial
#' and a no-go trial are the same statement about the same process, one
#' with the clock stopped where the experimenter stopped it and one
#' where the design stopped it.
#'
#' It goes through the variability average unconditionally, which is
#' exact rather than approximate when no variability is estimated: the
#' node sets are then the one-node rules, the widths are zero, and the
#' average is one evaluation at the mean, so this returns the same bits
#' `gng_lpdf()` returns for a no-go row at the same time.
#'
#' @noRd
gng_lccdf <- function(q, dpars, nd) {
  gv <- function(nm) if (is.null(dpars[[nm]])) 0 else dpars[[nm]]
  ddm_nogo_lprob_var(q, dpars[["mu"]], dpars[["bs"]], dpars[["bias"]],
                     dpars[["ndt"]], gv("sv"), gv("sz"), gv("st"), nd)
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

#' Fill in everything the family could not know until it had the data.
#'
#' Two quantities, both properties of the GO responses: the
#' non-decision-time bound, and the margin by which the non-decision
#' time range is held below each row's own response time. A no-go row's
#' entry is a placeholder, so neither is allowed to read one.
#'
#' The family is rebuilt rather than patched, because the `st` link is
#' scaled onto twice the bound and a link cannot be given its scale
#' after the fact. `ddm_ndt_finalize()` then installs the `ndt` link and
#' raises the bound refusal, so that the sentence a user sees is the one
#' the shared helper has always written.
#'
#' @noRd
gng_finalize <- function(cfg, y, aterms) {
  go <- gng_go_rt(y, aterms)
  lo <- min(go)
  fam <- gng_family(cfg, ub = cfg$max_ndt %||% lo, delta = 1e-9 * lo)
  ddm_ndt_finalize(fam, go, cfg$max_ndt, "wiener_gng")
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
  g <- function(nm, default) {
    if (is.null(dpars[[nm]])) rep_len(default, n)
    else rep_len(as.numeric(dpars[[nm]]), n)
  }
  mu <- g("mu", 0)
  bs <- g("bs", 1)
  ndt <- g("ndt", 0)
  bias <- g("bias", 0.5)
  # The variability widths reach the draw rather than the acceptance
  # test, because the process draws its per-trial parameters before it
  # runs and the acceptance test is on the OUTCOME. wiener() has to
  # reweight the drift by the boundary probability it implies; here the
  # rejection on "reached the go boundary before the deadline" does that
  # job already, and does it for all three at once.
  sv <- g("sv", 0)
  sz <- g("sz", 0)
  st <- g("st", 0)
  out <- td
  todo <- which(go > 0.5)
  for (round in 1:200) {
    if (!length(todo)) break
    d <- ddm_simulate(length(todo), mu = mu[todo], bs = bs[todo],
                      ndt = ndt[todo], bias = bias[todo],
                      sv = sv[todo], sz = sz[todo], st = st[todo])
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
#' @param sv,sz,st Across-trial variability in the drift rate, the
#'   relative start point and the non-decision time, as [ddm_simulate()]
#'   takes them. Zero, the default, is the plain model.
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
                                bias = 0.5, deadline = 1.5,
                                sv = 0, sz = 0, st = 0) {
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
  d <- ddm_simulate(n, mu = mu, bs = bs, ndt = ndt, bias = bias,
                    sv = sv, sz = sz, st = st)
  responded <- d$upper == 1 & d$rt <= td
  data.frame(responded = as.integer(responded),
             rt = ifelse(responded, d$rt, td),
             deadline = td)
}
