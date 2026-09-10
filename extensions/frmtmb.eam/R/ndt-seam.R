# The second thing this package offers another extension package: the
# BOUND a non-decision time is measured against. R/extension-api.R
# carries the first, the Wiener density itself.
#
# WHY IT EXISTS. The density is only half of what a joint
# choice-and-response-time family needs. The other half is the bound,
# because the density is zero at and below `ndt` and a fit walks over
# that edge unless the constraint is structural. frmtmb.learn's
# `rlddm()` wrote its own scaled logit for exactly that reason and,
# writing its own, inherited the defect item 1.0a removed here: one
# bound taken over the WHOLE data set, so a subject whose true
# non-decision time is above the global fastest response cannot be
# expressed at all.
#
# WHAT IS AND IS NOT PROMISED.
#
# Promised: the bound. How it is derived from the response, the
# `ndt_group()` addition term and `max_ndt`; the refusals that
# derivation owes; the link that makes it structural; the per-row floor
# reaching a density as data under the reserved name `ndt_floor`; and
# the record a fitted family carries, which `ndt_time()` and the
# unread-grouping frame check both read.
#
# NOT promised: any of the four densities' internals; the slot WRAPPING
# this package does for its own families, since `ddm_ndt_install()`
# rewrites `lpdf`, `lcdf`, `lccdf`, `sim` and `post$mean_fn` and that
# assumes this package's own family layout; the across-trial range `st`
# and its doubled bound; the label hash the floor table is keyed on; and
# the starting value of anything but `ndt`.
#
# WHERE THE BOUNDARY IS, AND WHY IT MOVED ONCE. `ddm_ndt_install()` is
# not exported and should not be: it rewrites five slots BY NAME, and
# the first consumer's likelihood is not in any of them. But the seam
# first stopped one function short of that, at the bound alone, and
# documented the multiplication as three lines a consumer copies. The
# arithmetic is three lines; the REFUSAL that goes with it is eleven,
# and it is the half that stops a fraction being read as seconds. So
# `ndt_apply()` is exported too. It takes a list in and returns a
# number, which assumes nothing about family layout, and it is the whole
# of what a consumer's density has to change.

#' The bound a non-decision time is measured against
#'
#' Derives the upper bound [wiener()] and its siblings put on the
#' non-decision time, so that a family in another package gets the same
#' bound, the same `ndt_group()` behavior and the same refusals without
#' writing them again. Pair it with [ndt_bound_attach()], which puts the
#' result on a family object.
#'
#' The bound is one of two things.
#'
#' * With no `ndt_group()` in `aterms` it is ONE number: `max_ndt` when
#'   that is given, and the fastest response otherwise. `ndt` is then a
#'   TIME, the link carries the bound, and a consumer's density reads
#'   `dpars$ndt` in the response's own units.
#' * With `ndt_group()` it is one number PER GROUP, that group's own
#'   fastest response. `ndt` is then a FRACTION of the row's own bound,
#'   on a plain logit, and the density multiplies it back out by the
#'   per-row vector that arrives in the addition-term values under the
#'   reserved name `ndt_floor`.
#'
#' `max_ndt` and `ndt_group()` together are refused: they set the same
#' bound to different things.
#'
#' @section Why the group's own floor and not the data set's:
#' The information about one subject's non-decision time is that
#' subject's own fastest response. At 30 subjects by 400 trials with a
#' between-subject spread of 26 ms on a mean of 250 ms, 20 of the 30
#' subjects have a true non-decision time above the GLOBAL minimum while
#' none is above its own, so a single bound runs the fit to the wall.
#' [wiener()] and `dev/ndt-findings.md` carry the measurements.
#'
#' @param y The response, a numeric vector of times.
#' @param aterms The addition-term values for this response, as a
#'   density receives them. Only `ndt_group` is read.
#' @param max_ndt An absolute bound in the response's own units, or
#'   `NULL` to take the fastest response. Above the fastest response it
#'   is refused, and it cannot be combined with `ndt_group()`.
#' @param what The family's name, used in the refusals so that the
#'   sentence a user reads names the family they wrote.
#'
#' @return An object of class `"frmtmb_eam_ndt_bound"`: a list with
#'   `ub`, the scalar bound; `floors`, the per-group bounds or `NULL`;
#'   `sizes`, the trial count behind each one; `pending`, always `FALSE`
#'   here; and `what`.
#' @seealso [ndt_bound_attach()] to put it on a family,
#'   [ndt_bound_of()] to read one back, [ndt_time()] to report a fitted
#'   non-decision time in the response's own units.
#' @examples
#' rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
#' g <- c("a", "a", "a", "b", "b")
#' # `[[` and not `$`: the record has five names and `bd$s` would
#' # partial-match `sizes`
#' ndt_bound(rt)[["ub"]]
#' ndt_bound(rt, list(ndt_group = ndt_bound_key(g)))[["floors"]]
#' @export
ndt_bound <- function(y, aterms = list(), max_ndt = NULL,
                      what = "ndt_bound()") {
  # `?ndt_bound` invites a caller to assemble `aterms` by hand and call
  # this outside frm(), and on that path there is no valid_y() in front
  # of it. Without these two checks `ndt_bound("a")` came back with
  # ub "a", and `ndt_bound(c(0.3, NA))` came back with ub NA_real_,
  # which is the exact sentinel ddm_scaled_logit() reads as PENDING: a
  # response with one missing value would then produce "the bound is not
  # set yet" from a family whose bound WAS set.
  ddm_check_y_times(y, what)
  ddm_check_max_ndt(max_ndt, what)
  sp <- ddm_ndt_spec(y, aterms, max_ndt, what)
  if (!is.null(max_ndt) && sp[["ub"]] > min(y)) {
    stop(what, ": max_ndt = ", format(sp[["ub"]]), " is above the ",
         "fastest response (", format(min(y)), "), so the bound admits ",
         "non-decision times at which that trial has no likelihood. ",
         "Give a max_ndt at or below the fastest response, or drop it ",
         "and let the data set the bound.", call. = FALSE)
  }
  ddm_ndt_bound_new(sp[["ub"]], sp[["floors"]], FALSE, sp[["sizes"]],
                    what)
}

#' The bound a family carries before it has seen a response
#'
#' A family object built outside [frmtmb::frm()] has no data, so it has
#' no fastest response and no bound. That state is not neutral, and this
#' is the seam entry that settles it, so that a family which reaches
#' [ndt_bound_attach()] at CONSTRUCTION is in one of two honest states
#' rather than carrying whatever link its constructor happened to
#' declare.
#'
#' * `max_ndt = NULL` gives a PENDING bound. The `ndt` link then refuses
#'   by name until [frmtmb::frm()] settles it. That refusal is
#'   load-bearing rather than decorative: [frmtmb::mixture()] never
#'   finalizes its components, so a component with an unsettled bound
#'   would read `ndt` on a scale nothing had set; and `bf(ndt = )` runs
#'   `linkfun()` at PARSE time, where the range check is against
#'   whichever link the family is carrying, so a constant outside the
#'   SETTLED bound passes a check against a placeholder and reaches the
#'   objective as `NaN`.
#' * `max_ndt = <a number>` gives a SETTLED bound, and the family scores
#'   outside `frm()`: `ndt` is a time measured against that bound,
#'   `predict(dpar = "ndt", type = "response")` reports seconds, and a
#'   pinned constant means what it says.
#'
#' The second is what a package reaching into another family's density
#' slot needs, and it is what [wiener()], [lba()], [rdm()] and
#' [wiener_gng()] do for themselves through their own `max_ndt`
#' argument.
#'
#' @param max_ndt An absolute bound in the response's own units, or
#'   `NULL` for the pending state.
#' @param what The family's name, for the refusals.
#' @return An object of class `"frmtmb_eam_ndt_bound"`.
#' @seealso [ndt_bound()] for the bound derived once the response is in
#'   hand, [ndt_bound_attach()] to put either on a family.
#' @examples
#' ndt_bound_pending()[["pending"]]
#' ndt_bound_pending(0.4)[["ub"]]
#' @export
ndt_bound_pending <- function(max_ndt = NULL, what = "ndt_bound()") {
  if (is.null(max_ndt)) {
    return(ddm_ndt_bound_new(NA_real_, NULL, TRUE, NULL, what))
  }
  ddm_check_max_ndt(max_ndt, what)
  ddm_ndt_bound_new(as.numeric(max_ndt), NULL, FALSE, NULL, what)
}

#' The one check both seam entry points make on `max_ndt`.
#'
#' Written once because the two disagreed: `ndt_bound_pending()` refused
#' a non-positive or non-scalar `max_ndt` by name and `ndt_bound()`
#' accepted `0`, `-1` and `NA`, which is two contracts on one argument
#' with one of them exported as API.
#'
#' @noRd
ddm_check_max_ndt <- function(max_ndt, what) {
  if (is.null(max_ndt)) return(invisible(NULL))
  if (!is.numeric(max_ndt) || length(max_ndt) != 1L ||
      !is.finite(max_ndt) || max_ndt <= 0) {
    stop(what, ": max_ndt is one positive number, the largest ",
         "non-decision time the fit may consider, in the response's ",
         "own units.", call. = FALSE)
  }
  invisible(NULL)
}

#' The response a bound can be derived from.
#'
#' The bound is the fastest response, so a response that is empty, not
#' numeric, not finite or not positive has no bound to give. `min()`
#' would otherwise return `Inf` on an empty vector, a character maximum
#' on a character vector, or `NA_real_` on one missing value, and the
#' last of those collides with the PENDING sentinel.
#'
#' @noRd
ddm_check_y_times <- function(y, what) {
  if (!is.numeric(y) || !length(y) || anyNA(y) || !all(is.finite(y)) ||
      any(y <= 0)) {
    stop(what, ": the response a non-decision-time bound is derived ",
         "from is one or more positive, finite times, and the bound is ",
         "the fastest of them. Got ", class(y)[[1L]], " of length ",
         length(y),
         if (is.numeric(y) && length(y) && anyNA(y)) {
           " with a missing value"
         } else if (is.numeric(y) && length(y) && any(y <= 0)) {
           " with a value at or below zero"
         } else "",
         ".", call. = FALSE)
  }
  invisible(NULL)
}

#' The non-decision time on the response's own scale
#'
#' The whole of what a per-group bound costs a consumer's density, and
#' it is one multiplication and one refusal. Call it wherever the
#' density, the mean or the simulator would have read `dpars$ndt`.
#'
#' With a per-group bound `dpars$ndt` is a FRACTION and the row's own
#' bound arrives in the addition-term values as `ndt_floor`, so this
#' returns `ndt * ndt_floor`. With a scalar bound `ndt_floor` is absent,
#' `dpars$ndt` is already a time, and this returns it untouched, so a
#' model that did not opt in is the arithmetic it always was.
#'
#' @section Why this is exported rather than documented:
#' The arithmetic is three lines. The REFUSAL is not, and it is the part
#' that matters: a grouped model whose density is reached without
#' `ndt_floor` would read the fraction as seconds and report a converged
#' fit at a non-decision time several times too small. That is the
#' silent wrong answer the per-group bound exists to remove, and a
#' consumer who copies the arithmetic and not the refusal reintroduces
#' it. The grouping itself is in the addition-term values on every path
#' that can reach a likelihood, so `ndt_group` present with `ndt_floor`
#' absent is exactly that condition, and this refuses on it.
#'
#' @section What happens if a consumer does not call this:
#' Nothing refuses, and the fit is silently wrong. This package cannot
#' check another package's arithmetic: a family that attaches a
#' per-group bound through [ndt_bound_attach()] and then reads
#' `dpars[["ndt"]]` directly is reading a fraction as a time, and
#' everything downstream still looks right.
#'
#' Measured, on two families identical in every line but this one, a
#' shifted gamma over four groups of 300 rows with true shifts 0.60 to
#' 0.90: both fit, both report convergence code 0, a positive definite
#' Hessian and no non-finite standard errors, both give the IDENTICAL
#' log-likelihood 1345.58425, and the one that skipped this call has an
#' [ndt_time()] up to **37.5 percent too small**. Every fitted time is
#' still below its own group's floor on both, which is why that check
#' cannot catch it: it tests the floor lookup, not the arithmetic.
#'
#' The fraction one half that [ndt_bound_attach()] sets as the `ndt`
#' starting value is a partial guard and only by luck. On a faster
#' design the same broken family died at once with `NA/NaN gradient
#' evaluation`, because half a second is past most of the data; on a
#' slower one it converged silently. A guard that works only when the
#' responses are fast is not a guard.
#'
#' @section What it does NOT do:
#' It returns a value rather than a modified `dpars`, because where the
#' value goes is the consumer's own layout. Two consequences a consumer
#' owns. A family with an across-trial RANGE on the non-decision time
#' scales that separately; this package's own families use twice the
#' bound. And a family that reads a linear predictor back through
#' frmtmb's robust dpar accessors must drop `dpars$.eta_ndt` itself
#' after calling this, because the linear predictor beside `ndt` no
#' longer maps to the value that was used.
#'
#' @param dpars The distributional parameters, as the density receives
#'   them. Only `ndt` is read.
#' @param aterms The addition-term values. Only `ndt_floor` and
#'   `ndt_group` are read. Defaults to `dpars`, for a family whose
#'   engine merges the two into one list, which is what
#'   `frmtmb.learn::rlddm()` does.
#' @param what The family's name, so that the refusal names the family
#'   the user wrote.
#'
#' @return The non-decision time in the response's own units, recycled
#'   to the length the arguments imply. An 'RTMB' advector passes
#'   through: there is no comparison and no branch on a parameter here,
#'   only two `is.null()` tests on data.
#' @seealso [ndt_bound_attach()], which puts `ndt_floor` where this
#'   reads it.
#' @examples
#' # no grouping: `ndt` is already a time
#' ndt_apply(list(ndt = 0.2))
#' # grouped: a fraction of each row's own bound
#' ndt_apply(list(ndt = 0.5), list(ndt_floor = c(0.4, 0.6)))
#' @export
ndt_apply <- function(dpars, aterms = dpars, what = "this family") {
  fl <- aterms[["ndt_floor"]]
  if (!is.null(fl)) return(dpars[["ndt"]] * fl)
  if (!is.null(aterms[["ndt_group"]])) {
    stop(what, ": this model bounds the non-decision time by each ",
         "ndt_group()'s own fastest response, so `ndt` is a fraction ",
         "of that bound and the per-row bound has to arrive with the ",
         "addition-term values as `ndt_floor`. It did not, and reading ",
         "the fraction as a time would give a converged fit at a wrong ",
         "answer. ndt_bound_attach() is what adds `ndt_floor`, so a ",
         "family object assembled by hand, or a frame built without ",
         "that family_finalize(), is what reaches this.", call. = FALSE)
  }
  dpars[["ndt"]]
}

#' The code an `ndt_group()` column is keyed on
#'
#' The addition-term registry coerces every `ndt_group()` column to a
#' numeric code before any family sees it, because an addition term's
#' value is baked into the tape as data. This is that coercion, exported
#' so that a caller holding a raw grouping can build the `aterms` list
#' [ndt_bound()] reads without going through [frmtmb::frm()].
#'
#' The code is a function of the group's LABEL and of nothing else, so a
#' factor, a character vector, a logical and integer codes that name the
#' same groups all give the same bound, and none of them depends on the
#' level order of the frame the column was evaluated in.
#'
#' @section It is a hash, and it is not injective:
#' The code is two polynomial rolling hashes over the label's code
#' points, packed into one double below `2^49`. Two different labels
#' CAN therefore share a code. Where both are in one column this refuses
#' by name rather than merging them. Where they are not, that is, where
#' a label the fit never saw collides with one it did, the unseen row is
#' paired with the collided group's bound instead of being refused, at a
#' rate of about 1.8e-15 per pairing.
#'
#' That residual is structurally zero for the identifier shapes anyone
#' uses: labels of one to three printable ASCII characters cannot
#' collide at all, and all 10^8 identifiers of the form `S00000000` hash
#' to 10^8 distinct codes where a birthday count predicts 8.88
#' collisions. Constructing an actual collision needed labels differing
#' in four positions over a span of 6401 code points.
#'
#' @param x A grouping: a factor, character, logical or numeric vector.
#' @return A numeric vector of codes, one per element of `x`.
#' @seealso [ndt_bound()], which reads the codes under the name
#'   `ndt_group`.
#' @examples
#' ndt_bound_key(factor(c("s1", "s2", "s1")))
#' @export
ndt_bound_key <- function(x) ddm_coerce_ndt_group(x)

#' Put a non-decision-time bound on a family
#'
#' Records `bound` on `fam`, sets the `ndt` link that makes the bound
#' structural, and, for a per-group bound, arranges for the per-row
#' floor to reach the family's likelihood as data. It does NOT touch the
#' likelihood: what a consumer's density has to do is three lines,
#' written below, and doing it there rather than through a wrapper is
#' what keeps this seam free of every assumption about where that
#' density lives.
#'
#' @section What the consumer's density must do:
#' Call [ndt_apply()] wherever it would have read `dpars[["ndt"]]`, in
#' the density and in the mean and in the simulator alike. A family that
#' attaches a bound here and then does NOT call it fits, converges and
#' reports a non-decision time measured at 37.5 percent too small, with
#' nothing refusing; `?ndt_apply` has that measurement. Nothing on this
#' side can check it. That is the
#' whole of it, and it is one line rather than three, because the three
#' lines of arithmetic come with a refusal that is longer than they are
#' and that is the part that matters. An earlier version of this seam
#' documented the arithmetic and left the refusal in prose; a consumer
#' who copied the first and not the second got a fraction silently read
#' as seconds, which is the defect the per-group bound exists to remove.
#'
#' @section What it changes on the family:
#' `links$ndt`, `ndt_bound`, and for a per-group bound `aterm_data` and
#' the `ndt` starting value, which becomes the fraction one half. A
#' scalar bound leaves the family's own starting value alone. It also
#' records the family's pre-attach `aterm_data` and `ndt` starting value
#' in a reserved slot named `ndt_seam`, which is what makes a second
#' call replace the first bound instead of composing with it; a family
#' with a slot of that name of its own would lose it. Nothing
#' else moves, and `st` is not handled: a family with an across-trial
#' range on the non-decision time sets that bound itself.
#'
#' Calling it twice replaces the first bound rather than stacking on it,
#' which matters because [frmtmb::frm()] runs `family_finalize` again on
#' the `influence()`, `frm_simulate()` and prior-predictive paths.
#'
#' @param fam A `frmtmb_family` with a distributional parameter named
#'   `ndt`.
#' @param bound The result of [ndt_bound()].
#' @return `fam`, modified.
#' @seealso [ndt_bound()], [ndt_bound_of()].
#' @examples
#' # what a family_finalize() written against this seam looks like
#' finalize <- function(fam, y, aterms) {
#'   if (!is.null(ndt_bound_of(fam))) return(fam)
#'   ndt_bound_attach(fam, ndt_bound(y, aterms, what = "my_family"))
#' }
#' @export
ndt_bound_attach <- function(fam, bound) {
  if (!inherits(fam, "frmtmb_family")) {
    stop("ndt_bound_attach(fam =) takes the family object a ",
         "family_finalize() is handed, and got ",
         paste(class(fam), collapse = "/"), call. = FALSE)
  }
  if (!inherits(bound, "frmtmb_eam_ndt_bound")) {
    stop("ndt_bound_attach(bound =) takes what ndt_bound() returns. A ",
         "hand-built list is refused because the seam owns the ",
         "refusals that go with deriving a bound, and a list assembled ",
         "beside them would carry none of them.", call. = FALSE)
  }
  # THE GUARD, AND WHY IT IS NOT `ndt_raw`. The first version refused a
  # family carrying `ndt_raw`, the marker ddm_ndt_install() leaves. Only
  # FOUR of this package's five families carry it. gddm() has an `ndt`,
  # has no `ndt_raw`, and its density reads no `ndt_floor` at all, so it
  # passed: the attach replaced its bounded scaled logit with a plain
  # logit, planted an `ndt_floor` nothing reads, and moved its starting
  # value from 0.155 s to 0.5. The guard tested a slot that happened to
  # correlate with the property instead of the property.
  #
  # `ddm_accepts` is the table this package's five families are
  # enumerated in, and the compatibility rows already derive from it, so
  # a sixth family added there is refused here without anyone
  # remembering. `ndt_raw` stays as a second condition, for a family
  # object that has already been through ddm_ndt_install().
  if (!is.null(fam[["ndt_raw"]]) ||
        isTRUE(fam[["family"]] %in% names(ddm_accepts))) {
    stop("ndt_bound_attach(): `", fam[["family"]], "` is one of ",
         "frmtmb.eam's own families and cannot take a bound this way. ",
         "wiener(), lba(), rdm() and wiener_gng() install theirs ",
         "through their constructor's max_ndt argument and their own ",
         "family_finalize(), and attaching a second one here would set ",
         "the link without the slot wrapping that goes with it, so the ",
         "density would read a fraction as a time. gddm() takes no ",
         "per-group bound at all: its solver reads every distributional ",
         "parameter at the FIRST ROW of a condition, so a per-row bound ",
         "would never reach its density, and ?gddm says so. This seam ",
         "is for a family in ANOTHER package whose own density reads ",
         "`ndt_floor`.", call. = FALSE)
  }
  if (!"ndt" %in% names(fam[["links"]])) {
    stop("ndt_bound_attach(): this family has no distributional ",
         "parameter named `ndt`, so there is nothing to bound. The ",
         "seam sets the `ndt` link and the `ndt` starting value and ",
         "reads no other parameter.", call. = FALSE)
  }
  # The pre-attach state is kept so that a second call REPLACES the
  # first bound instead of composing with it. frm() finalizes the same
  # family again on the influence(), frm_simulate() and prior-predictive
  # paths, and a second composed aterm_data() would put two `ndt_floor`
  # entries in the addition-term values.
  if (is.null(fam[["ndt_seam"]])) {
    fam[["ndt_seam"]] <- list(aterm_data = fam[["aterm_data"]],
                              init_ndt = fam[["init_dpars"]][["ndt"]])
  }
  prev <- fam[["ndt_seam"]][["aterm_data"]]
  fam[["aterm_data"]] <- prev
  fam[["init_dpars"]][["ndt"]] <- fam[["ndt_seam"]][["init_ndt"]]
  fl <- bound[["floors"]]
  what <- bound[["what"]]
  if (is.null(what)) what <- "ndt_bound()"
  if (is.null(fl)) {
    fam[["links"]][["ndt"]] <- ddm_scaled_logit(
      if (isTRUE(bound[["pending"]])) NA_real_ else bound[["ub"]],
      "non-decision time", what)
  } else {
    scale <- ddm_ndt_scaler(fl, what)
    fam[["links"]][["ndt"]] <- "logit"
    # half of each group's own bound, which is what the scalar case's
    # `0.5 * min(y)` meant when there was one bound. max_ndt cannot be
    # in play here, because the pair is refused, so half the bound is
    # always a time some row reached.
    fam[["init_dpars"]][["ndt"]] <- function(y, aterms) 0.5
    fam[["aterm_data"]] <- if (is.null(prev)) {
      function(y, aterms) list(ndt_floor = scale(aterms))
    } else {
      function(y, aterms) {
        c(prev(y, aterms), list(ndt_floor = scale(aterms)))
      }
    }
  }
  fam[["ndt_bound"]] <- bound
  fam
}

#' The non-decision-time bound a family already carries
#'
#' `NULL` when the family has none, and `NULL` when the one it has is
#' still PENDING, which is the state a family constructed outside
#' [frmtmb::frm()] is in before the response has been seen.
#'
#' A settled bound is KEPT rather than re-derived, and that is what this
#' is for. [frmtmb::frm()] runs `family_finalize` again on new rows for
#' `influence()`, `frm_simulate()` and the prior-predictive path, so a
#' bound re-derived from fewer rows would make a leave-one-out refit a
#' refit of a DIFFERENT model. A `family_finalize` written against this
#' seam therefore asks first and returns the family untouched when the
#' answer is not `NULL`.
#'
#' @param fam A `frmtmb_family` object, fitted or not.
#' @return An object of class `"frmtmb_eam_ndt_bound"`, or `NULL`.
#' @seealso [ndt_bound()], [ndt_bound_attach()].
#' @examples
#' # a family built outside frm() has no settled bound yet
#' ndt_bound_of(wiener())
#' ndt_bound_of(wiener(max_ndt = 0.4))[["ub"]]
#' @export
ndt_bound_of <- function(fam) {
  bd <- fam[["ndt_bound"]]
  if (is.null(bd) || isTRUE(bd[["pending"]])) return(NULL)
  bd
}
