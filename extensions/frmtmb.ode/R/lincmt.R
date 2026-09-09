# The analytic one- to three-compartment linear model.
#
# frm_ode() solves the system numerically, one adjoint solve per group
# per segment, and Phase 0 measured that at 3948 s on a 100 subject
# design (dev/scale-findings.md, the `ode` row). A LINEAR compartment
# model does not need a solver: the state is a superposition of the
# impulse responses of the doses that came before it, and each impulse
# response is a sum of exponentials whose rate constants are the
# eigenvalues of the rate matrix. frm_lincmt() writes that sum.
#
# The whole difficulty is that the textbook form of such a sum divides
# by differences of rate constants, and those differences go to zero on
# a set an optimizer walks through: ka equal to ke is the ordinary
# flip-flop of an oral model, and it is the starting value of any fit
# whose two log rates start at the same number. Every primitive below is
# therefore written so that no difference of rate constants is ever
# formed in a denominator, and so that no intermediate overflows. See
# dev/lincmt-findings.md for the measurements that back that up.
#
# The schedule grammar is frm_ode()'s, and it is not re-implemented: the
# `events` table goes through ode_split_events(), so `ii`, `addl`, `ss`,
# `duration`, `method` and `group` mean exactly what ?frm_ode says they
# mean. What frm_lincmt() cannot express it refuses by name.

# A zero denominator cannot be branched around on the automatic
# differentiation tape, because RTMB refuses comparison on AD types. It
# is offset instead. The offset is below the smallest rate constant any
# real arithmetic produces, so it changes no answer and it turns a 0/0
# into a finite number. It is 1e-150 rather than something smaller so
# that its SQUARE, which a first derivative of `x / z` puts in a
# denominator, is still a normal double.
lincmt_tiny <- 1e-150

#' `max(x, 0)`, without a comparison.
#'
#' @noRd
lincmt_relu <- function(x) (x + abs(x)) / 2

# Taylor coefficients of (1 - exp(-x)) / x, that is (-1)^k / (k + 1)!,
# to the order that leaves 1e-19 at x = 0.5.
lincmt_phi_coef <- (-1)^(0:12) / factorial(1:13)

#' `(1 - exp(-x)) / x`, for `x >= 0`, equal to 1 at zero.
#'
#' The offset alone is not enough, and finding out why is the substance
#' of this file. `-expm1(-z) / z` at `z = 1e-150` VALUES correctly, but
#' its derivative is taken by the chain rule as `exp(-z)/z - n/z^2`,
#' two numbers of size `1e150` whose difference is `-1/2`: every digit
#' of the gradient is lost, and the tape returned 1.2e286 where the
#' answer is 6.04. A model started with two log rates at the same value
#' hits that on its first evaluation.
#'
#' So the series is used where it converges and the `expm1()` form
#' where the series would need too many terms, blended by a smoothstep
#' whose first and second derivatives vanish at both ends. Nothing here
#' branches: both arms are evaluated for every element, the series
#' argument is capped so that it cannot overflow in the arm whose
#' weight is zero, and the weight is built from `abs()`.
#'
#' @noRd
lincmt_phi <- function(x) {
  z <- x + lincmt_tiny
  naive <- -expm1(-z) / z
  # capped so that x^12 cannot overflow where the weight is zero anyway
  xs <- x - lincmt_relu(x - 1)
  s <- lincmt_phi_coef[[13L]]
  for (k in 12:1) s <- lincmt_phi_coef[[k]] + xs * s
  t <- lincmt_relu(x - 0.125) / 0.375
  t <- t - lincmt_relu(t - 1)
  w <- 1 - t * t * t * (10 + t * (-15 + 6 * t))
  w * s + (1 - w) * naive
}

#' `(exp(-p) - exp(-q)) / (q - p)` for `p, q >= 0`.
#'
#' The building block of every response below. It is symmetric in its
#' arguments and it never overflows: the exponential that is factored
#' out is the LARGER of the two, chosen with `abs()` rather than with a
#' comparison, so `exp()` is only ever called on a non-positive number.
#' The naive spelling loses all precision as `p` approaches `q` and
#' returns `NaN` when they are equal; this one returns `exp(-p)`.
#'
#' @noRd
lincmt_diff <- function(p, q) {
  d <- abs(p - q)
  exp(-(p + q - d) / 2) * lincmt_phi(d)
}

#' `E2(a, b; u)`: the convolution of two exponentials.
#'
#' `integral over [0, u] of exp(-a s) exp(-b (u - s)) ds`, which is
#' `(exp(-a u) - exp(-b u)) / (b - a)` and `u exp(-a u)` when `a == b`.
#' This is the whole first-order absorption term.
#'
#' @noRd
lincmt_e2 <- function(a, b, u) u * lincmt_diff(a * u, b * u)

#' `1 / (1 - exp(-lam * ii))`: the accumulation factor of a dose
#' repeated forever every `ii`.
#'
#' @noRd
lincmt_geo <- function(lam, ii) -1 / expm1(-(lam * ii + lincmt_tiny))

#' The sum of `E2(a, b; u + j ii)` over every `j >= 0`.
#'
#' The steady state of a repeated dose through a first-order absorption
#' step. Summing the textbook partial fractions term by term would
#' divide by `b - a` twice; putting the two geometric series over a
#' common denominator first leaves a numerator that is two `lincmt_diff`
#' calls, so the singularity at `a == b` is removed rather than
#' approached.
#'
#' @noRd
lincmt_e2_ss <- function(a, b, u, ii) {
  num <- u * lincmt_diff(a * u, b * u) +
    (ii - u) * lincmt_diff(a * u + b * ii, b * u + a * ii)
  num * lincmt_geo(a, ii) * lincmt_geo(b, ii)
}

#' Eigenvalues of the disposition matrix, and the coefficients of the
#' central compartment's impulse response.
#'
#' Returns `lam` and `coef` such that the amount left in the central
#' compartment `u` after a unit bolus into it is
#' `sum(coef[[i]] * exp(-lam[[i]] * u))`. For a mammillary model the
#' coefficients are non-negative and sum to one, so the response is a
#' positive combination and nothing in it cancels; that is why the
#' textbook cancellation shows up in the ABSORPTION term and not here.
#'
#' Two places still need care and get it. The smaller root of the
#' two-compartment quadratic is taken from the PRODUCT of the roots,
#' because `m - Delta` cancels when the peripheral compartment is nearly
#' decoupled. The cubic is solved trigonometrically, with its two
#' arguments clamped by `abs()` rather than by a comparison, because
#' rounding can push a discriminant that is zero in exact arithmetic
#' just outside the domain of `sqrt()` or `acos()`.
#'
#' Every argument may be a vector, one element per group.
#'
#' @noRd
lincmt_disp <- function(ncmt, ke, k12, k21, k13, k31) {
  if (ncmt == 1L) return(list(lam = list(ke), coef = list(1 + 0 * ke)))
  if (ncmt == 2L) {
    g <- (ke + k12 - k21) / 2
    dd <- sqrt(g * g + k12 * k21 + lincmt_tiny)
    lam1 <- (ke + k12 + k21) / 2 + dd
    lam2 <- ke * k21 / lam1
    c1 <- (dd + g) / (2 * dd)
    return(list(lam = list(lam1, lam2), coef = list(c1, 1 - c1)))
  }
  a2 <- ke + k12 + k13 + k21 + k31
  a1 <- ke * k21 + ke * k31 + k21 * k31 + k12 * k31 + k13 * k21
  a0 <- ke * k21 * k31
  # Every one of `ke`, `k12`, `k21`, `k13` and `k31` underflowing to
  # exactly zero makes `pp` and `qq` both zero and returns NaN from
  # `3 qq / (2 pp)`. That is left as a NaN deliberately: offsetting
  # `pp` removes the NaN and returns a finite WRONG answer, because the
  # eigenvalues are then of order 1e-75 and `lincmt_tiny` is no longer
  # negligible against their differences, so the coefficients come out
  # 4/9 and 4/9 instead of 1/2 and 1/2 and the whole trajectory is
  # eight-ninths of its true height. A silent factor of 1.125 is worse
  # than a NaN, and reaching this needs five log rates below -745 at
  # once. Measured by the review in dev/rev-lincmt-n1.R; this lane's
  # own dev/lincmt/lincmt-n1.R shows the NaN and that ordinary rates
  # are unaffected, but does NOT apply the offset.
  pp <- a1 - a2 * a2 / 3
  qq <- -2 * a2^3 / 27 + a1 * a2 / 3 - a0
  rr <- sqrt(lincmt_relu(-pp / 3) + lincmt_tiny)
  arg <- (3 * qq) / (2 * pp) * sqrt(lincmt_relu(-3 / pp))
  arg <- arg / sqrt(1 + lincmt_relu(arg * arg - 1))
  th <- acos(arg) / 3
  s <- a2 / 3
  l1 <- 2 * rr * cos(th) + s
  l2 <- 2 * rr * cos(th - 2 * pi / 3) + s
  l3 <- 2 * rr * cos(th - 4 * pi / 3) + s
  cf <- function(x, y, z) {
    (x - k21) * (x - k31) / ((x - y) * (x - z) + lincmt_tiny)
  }
  list(lam = list(l1, l2, l3),
       coef = list(cf(l1, l2, l3), cf(l2, l1, l3), cf(l3, l1, l2)))
}

# ---------------------------------------------------------------------
# Parameters.

# The two spellings a pharmacokinetic model is written in. Which one is
# meant is read off the names present, and a table that mixes them is
# refused rather than half-read.
lincmt_k_names <- c("ke", "V", "ka", "k12", "k21", "k13", "k31")
lincmt_cl_names <- c("CL", "V", "ka", "Q2", "V2", "Q3", "V3")
# `V` and `ka` are spelled the same in both, so they decide nothing
lincmt_shared_names <- c("V", "ka")

#' Which parameters this model shape needs, in each spelling.
#'
#' @noRd
lincmt_needed <- function(ncmt, depot, kind, want_v) {
  if (kind == "k") {
    nm <- if (want_v) c("ke", "V") else "ke"
    if (ncmt >= 2L) nm <- c(nm, "k12", "k21")
    if (ncmt == 3L) nm <- c(nm, "k13", "k31")
  } else {
    nm <- c("CL", "V")
    if (ncmt >= 2L) nm <- c(nm, "Q2", "V2")
    if (ncmt == 3L) nm <- c(nm, "Q3", "V3")
  }
  if (depot) nm <- c(nm, "ka")
  nm
}

#' One column, named, or an error naming the argument.
#'
#' @noRd
lincmt_col <- function(x, n_obs, arg, nm) {
  cols <- ode_columns(x, n_obs, arg)
  if (length(cols) != 1L) {
    stop("`", arg, "$", nm, "` has ", length(cols), " columns; each ",
         "element is ONE value per observation, constant within group, ",
         "or one value shared by every group", call. = FALSE)
  }
  cols[[1L]]
}

#' Validate `parms` and return the micro rate constants.
#'
#' @noRd
lincmt_rates <- function(pl, ncmt, depot, want_v) {
  if (!is.list(pl) || inherits(pl, "advector") || !length(pl)) {
    stop("`parms` must be a named list, for example ",
         "parms = list(ka = exp(lka), ke = exp(lke), V = exp(lV))",
         call. = FALSE)
  }
  nm <- names(pl)
  if (is.null(nm) || any(!nzchar(nm))) {
    stop("every element of `parms` must be named. frm_lincmt() reads ",
         "its parameters by name, not by position, because which ones ",
         "a model needs depends on `ncmt` and `depot`. The names are: ",
         paste(lincmt_k_names, collapse = ", "), " (rate constants), ",
         "or ", paste(lincmt_cl_names, collapse = ", "),
         " (clearances)", call. = FALSE)
  }
  if (anyDuplicated(nm)) {
    stop("`parms` names ", nm[anyDuplicated(nm)], " more than once",
         call. = FALSE)
  }
  has_k <- any(nm %in% setdiff(lincmt_k_names, lincmt_shared_names))
  has_cl <- any(nm %in% setdiff(lincmt_cl_names, lincmt_shared_names))
  if (has_k && has_cl) {
    stop("`parms` mixes the two parameterizations: ",
         paste(setdiff(intersect(nm, lincmt_k_names),
                       lincmt_shared_names), collapse = ", "),
         " are rate constants and ",
         paste(setdiff(intersect(nm, lincmt_cl_names),
                       lincmt_shared_names), collapse = ", "),
         " are clearances and volumes. Write one or the other; a mixed ",
         "table has no reading that is not a guess", call. = FALSE)
  }
  kind <- if (has_cl) "cl" else "k"
  want <- lincmt_needed(ncmt, depot, kind, want_v || kind == "cl")
  miss <- setdiff(want, nm)
  if (length(miss)) {
    stop("`parms` is missing ", paste(miss, collapse = ", "),
         ". A ", ncmt, "-compartment model ",
         if (depot) "with" else "without", " a depot needs ",
         paste(want, collapse = ", "), call. = FALSE)
  }
  # `V` is always allowed, so that changing `output` does not change
  # what `parms` must hold; it is only REQUIRED where it is read
  extra <- setdiff(nm, c(want, "V"))
  if (length(extra)) {
    stop("`parms` has ", paste(extra, collapse = ", "),
         ", which a ", ncmt, "-compartment model ",
         if (depot) "with" else "without",
         " a depot does not use. It needs exactly ",
         paste(want, collapse = ", "),
         ". Raise `ncmt`, or set depot = TRUE, or drop the parameter",
         call. = FALSE)
  }
  pl
}

# ---------------------------------------------------------------------
# The response of one compartment to one term of the schedule.

#' The amount in the output compartment, one value per term.
#'
#' `kind` is 1 for a single dose and 2 for a dose repeated every `ii`
#' back to minus infinity; `cmt` is 1 for the depot and 2 for the
#' central compartment. Every one of those is data, so the branches
#' below are taken once, when the tape is built, and each branch is one
#' vectorized expression over the terms that took it.
#'
#' @noRd
lincmt_resp <- function(u, dur, ii, kind, cmt, lam, coef, ka, out) {
  "[<-" <- RTMB::ADoverload("[<-")
  v <- numeric(length(u))
  nc <- length(lam)
  at <- function(x, i) if (length(x) == 1L) x else x[i]
  # sum over the disposition exponentials, for one subset of terms
  over <- function(i, f) {
    s <- f(1L, i)
    if (nc > 1L) for (k in 2L:nc) s <- s + f(k, i)
    s
  }
  if (out == 1L) {
    # the depot is fed by nothing, so only a depot dose reaches it
    i <- which(cmt == 1L & kind == 1L)
    if (length(i)) v[i] <- exp(-at(ka, i) * u[i])
    i <- which(cmt == 1L & kind == 2L)
    if (length(i)) {
      v[i] <- exp(-at(ka, i) * u[i]) * lincmt_geo(at(ka, i), ii[i])
    }
    return(v)
  }
  i <- which(cmt == 1L & kind == 1L)
  if (length(i)) {
    v[i] <- at(ka, i) * over(i, function(k, i)
      at(coef[[k]], i) * lincmt_e2(at(ka, i), at(lam[[k]], i), u[i]))
  }
  i <- which(cmt == 1L & kind == 2L)
  if (length(i)) {
    v[i] <- at(ka, i) * over(i, function(k, i)
      at(coef[[k]], i) *
        lincmt_e2_ss(at(ka, i), at(lam[[k]], i), u[i], ii[i]))
  }
  i <- which(cmt == 2L & kind == 1L & dur == 0)
  if (length(i)) {
    v[i] <- over(i, function(k, i)
      at(coef[[k]], i) * exp(-at(lam[[k]], i) * u[i]))
  }
  # An infusion delivers `value` over `duration`, so its response per
  # unit AMOUNT carries the 1 / duration of the rate. Still running at
  # the observation:
  i <- which(cmt == 2L & kind == 1L & dur > 0 & u <= dur)
  if (length(i)) {
    v[i] <- over(i, function(k, i)
      at(coef[[k]], i) * (u[i] / dur[i]) *
        lincmt_phi(at(lam[[k]], i) * u[i]))
  }
  # and finished:
  i <- which(cmt == 2L & kind == 1L & dur > 0 & u > dur)
  if (length(i)) {
    v[i] <- over(i, function(k, i)
      at(coef[[k]], i) * lincmt_phi(at(lam[[k]], i) * dur[i]) *
        exp(-at(lam[[k]], i) * (u[i] - dur[i])))
  }
  i <- which(cmt == 2L & kind == 2L & dur == 0)
  if (length(i)) {
    v[i] <- over(i, function(k, i)
      at(coef[[k]], i) * exp(-at(lam[[k]], i) * u[i]) *
        lincmt_geo(at(lam[[k]], i), ii[i]))
  }
  # a repeated infusion: `duration <= ii` is checked by the schedule, and
  # the nearest repeat is at least `ii` back, so every term is past its
  # own infusion
  i <- which(cmt == 2L & kind == 2L & dur > 0)
  if (length(i)) {
    v[i] <- over(i, function(k, i)
      at(coef[[k]], i) * lincmt_phi(at(lam[[k]], i) * dur[i]) *
        exp(-at(lam[[k]], i) * (u[i] - dur[i])) *
        lincmt_geo(at(lam[[k]], i), ii[i]))
  }
  v
}

# ---------------------------------------------------------------------
# What the closed form cannot express.

#' Refuse the parts of the `events` grammar superposition cannot carry.
#'
#' Each of these is a model frm_ode() fits and frm_lincmt() does not, so
#' each is named rather than approximated.
#'
#' @noRd
lincmt_check_events <- function(ev, label, n_dep, n_state) {
  bad <- setdiff(unique(ev[["method"]]), c("add", "reset"))
  if (length(bad)) {
    stop("`events` has ", paste(bad, collapse = " and "),
         " rows, which frm_lincmt() does not fit. The closed form is a ",
         "superposition of the doses that came before an observation, ",
         "and \"replace\" and \"multiply\" are not doses: they need the ",
         "state of every compartment at the instant they act, including ",
         "the peripheral compartments, which superposition never forms. ",
         "Use frm_ode()", call. = FALSE)
  }
  rst <- ev[["method"]] == "reset"
  if (any(rst & ev[["value"]] != 0)) {
    stop("`events` has a \"reset\" row with a value of ",
         format(ev[["value"]][rst & ev[["value"]] != 0][[1L]]),
         ". frm_lincmt() reads a reset to ZERO, which is NONMEM's and ",
         "rxode2's EVID = 3, as \"forget every dose before this time\", ",
         "and that needs no state. A reset to a non-zero level sets the ",
         "peripheral compartments as well. Use frm_ode()", call. = FALSE)
  }
  cm <- ev[["state"]][!rst]
  if (length(cm) && any(cm > n_dep + 1L)) {
    stop("`events` doses compartment ", max(cm), " of group '", label,
         "', which is a peripheral compartment. frm_lincmt() doses the ",
         if (n_dep) "depot or the central compartment" else
           "central compartment",
         " only: a peripheral compartment has no route of ",
         "administration, and its impulse response is the one piece of ",
         "the closed form that has no cancellation-free spelling. Use ",
         "frm_ode()", call. = FALSE)
  }
  if (n_dep && any(ev[["duration"]] > 0 & ev[["state"]] == 1L)) {
    stop("`events` infuses into the depot of group '", label,
         "'. That is zero-order absorption, and frm_lincmt() does not ",
         "fit it: the three-node convolution it needs cannot be written ",
         "without dividing by a difference of rate constants. An ",
         "infusion into the CENTRAL compartment is supported. Use ",
         "frm_ode()", call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------
# The exported helper.

#' Analytic one- to three-compartment pharmacokinetics
#'
#' Evaluates a linear compartment model in closed form, as an
#' alternative to [frm_ode()] for the systems where one exists. The
#' dosing grammar, the grouping and the place in a `bf(nl = TRUE)`
#' formula are [frm_ode()]'s, so a linear model written for the solver
#' becomes a `frm_lincmt()` call with the same `events` table. There is
#' no integrator: the states are a superposition of the impulse
#' responses of the doses that precede each observation, and a
#' steady-state record is the limit of that superposition rather than a
#' run-in.
#'
#' @details
#' # The model
#'
#' One central compartment, up to two peripheral compartments exchanging
#' with it, and an optional depot feeding it by first-order absorption:
#'
#' \deqn{
#'   \frac{dA_d}{dt} = -k_a A_d, \qquad
#'   \frac{dA_1}{dt} = k_a A_d - (k_e + k_{12} + k_{13}) A_1
#'                     + k_{21} A_2 + k_{31} A_3
#' }
#' \deqn{
#'   \frac{dA_j}{dt} = k_{1j} A_1 - k_{j1} A_j
#' }
#'
#' The compartments are named `"depot"`, `"central"`, `"peripheral1"`
#' and `"peripheral2"`, in that order, and `events$state` may use those
#' names or the positions they stand at. States are AMOUNTS; the default
#' `output = "conc"` divides the central amount by `V`.
#'
#' # Parameters
#'
#' `parms` is a NAMED list, because which parameters a model needs
#' depends on `ncmt` and `depot`. Two spellings are read, and mixing
#' them is refused:
#'
#' - rate constants: `ke`, `V`, and `ka`, `k12`, `k21`, `k13`, `k31` as
#'   the shape requires;
#' - clearances: `CL`, `V`, and `ka`, `Q2`, `V2`, `Q3`, `V3`, with
#'   \eqn{k_e = CL/V}, \eqn{k_{12} = Q_2/V}, \eqn{k_{21} = Q_2/V_2}.
#'
#' Each element is one value per observation, constant within a group,
#' or a single value shared by every group, exactly as `frm_ode()`'s
#' `parms` is. They are ordinary nonlinear parameters, so they take
#' fixed effects, random effects and covariates.
#'
#' # Dosing
#'
#' `events` is [frm_ode()]'s table and goes through the same validation,
#' so `time`, `value`, `state`, `method`, `duration`, `ii`, `addl`, `ss`
#' and `group` mean what `?frm_ode` says they mean, including that an
#' observation at a dose time reads the trough. `event_scale` is the
#' same estimated multiplier on every dose, which is how a
#' bioavailability is written.
#'
#' A steady-state row is exact here. `frm_ode()` approaches a steady
#' state by simulating `n_ss` cycles; `frm_lincmt()` sums the geometric
#' series, so the default `n_ss = Inf` is the limit itself and costs one
#' term. A finite `n_ss` writes out that many cycles as ordinary doses,
#' which is what reproduces a `frm_ode()` fit exactly and is otherwise
#' not worth paying for. The limit is what it says: a drug that is not
#' eliminated over its dosing interval accumulates without bound, and
#' the series says so, where a run-in of twenty cycles reports twenty
#' cycles worth. Nothing overflows, so an optimizer that wanders into
#' that region gets a large finite objective and backs out of it.
#'
#' # What is refused, and why
#'
#' Superposition adds impulse responses. It never forms the state
#' vector, so anything that reads or sets that vector is refused by
#' name rather than approximated:
#'
#' - `method = "replace"` and `method = "multiply"`. A `"reset"` to zero
#'   IS supported, because "forget every dose before this time" needs no
#'   state; a reset to a non-zero level is refused.
#' - a dose into a peripheral compartment, and `output` naming one.
#' - an infusion into the depot, which is zero-order absorption.
#' - `tv` and `tv_break`. A rate constant that changes with time makes
#'   the system time-varying, and superposition over the whole history
#'   is then wrong rather than approximate.
#' - an infusion that is still running at a `ss` or `"reset"` time.
#'
#' [frm_ode()] fits every one of them.
#'
#' # Accuracy
#'
#' The textbook closed form divides by differences of rate constants and
#' loses all its precision where two of them meet: \eqn{k_a = k_e} is
#' the ordinary flip-flop case of an oral model, and it is also the
#' starting value of any fit that starts two log rates at the same
#' number. **The ABSORPTION cancellation is removed completely**: every
#' response is built from `(exp(-p) - exp(-q)) / (q - p)` evaluated as
#' `exp(-max(p, q)) (1 - exp(-|p - q|)) / |p - q|`, which is exact at
#' `p == q` and never overflows, where the textbook difference of
#' exponentials returns `NaN` there and has lost most of its digits by
#' a separation of 1e-12. The disposition coefficients of a mammillary
#' model are non-negative and sum to one, so the sum of exponentials is
#' a positive combination.
#'
#' One cancellation is NOT removed, and it is worth knowing which. The
#' smaller of the two disposition coefficients of a two-compartment
#' model is spelled `1 - c1`, and when the peripheral compartment is
#' nearly decoupled that subtraction loses about as many digits as
#' `k12 k21 / g^2` has decades, much as the textbook spelling does. It
#' is bounded where it matters and unbounded where it does not, and the
#' two statements are the same statement: the ABSOLUTE error is held at
#' a few parts in 1e10 of the trajectory's own maximum, so the
#' POINTWISE relative error at a point `10^-d` below `Cmax` is at most
#' about that times `10^d`. An observation within four decades of
#' `Cmax` therefore carries at most about 1e-11 and the 1e-8 the
#' feature is held to is crossed only beyond seven decades below
#' `Cmax`, which is far below any assay's limit of quantification.
#'
#' Measured: over 118 schedules the disagreement with [frm_ode()] run
#' at `atol = rtol = 1e-12` is at most 2.6e-12 of the trajectory's own
#' scale, and tightening the solver drives it down rather than leaving
#' it, so it is the solver's tolerance rather than the closed form's
#' error. Against a 240-bit reference that never forms an eigenvalue,
#' the worst relative error through every coalescence down to exact
#' equality is 7.2e-15, and over a seven-decade box of rate constants
#' the worst error relative to the trajectory's maximum found by random
#' search is 3.6e-10. The repository's `dev/lincmt-findings.md` has the
#' tables.
#'
#' # Sampling
#'
#' `frmtmb.sample::frm_sample()` refuses a model containing
#' [frm_ode()], because \pkg{RTMBode} calls \pkg{deSolve} unguarded and
#' the chain aborts at warmup iteration 1. `frm_lincmt()` calls
#' neither, carries no refusal in [frmtmb::frm_compat()], and samples.
#'
#' # Cost
#'
#' One term per (observation, dose) pair, plus one term per group for a
#' steady-state record, all evaluated in one vectorized pass over the
#' whole data. There is no per-group loop on the tape and no solver, so
#' the cost is linear in the number of doses each observation follows.
#' A schedule with hundreds of `addl` doses before each observation is
#' the case where `frm_ode()`'s per-segment cost can win.
#'
#' @section The accurate domain:
#' **Up to a six-decade spread of rate constants, nothing degrades.**
#' Over 60 random three-compartment draws per cell, the worst
#' disagreement with [frm_ode()] at `atol = rtol = 1e-12` is 9.9e-11 at
#' a two-decade spread, 1.3e-09 at four and 1.4e-07 at six, and in
#' every one of those cells it is BELOW what the solver's own tolerance
#' costs on the same draws (1.3e-07, 1.6e-06, 1.9e-05). The direction
#' an optimizer actually travels is cleaner still: driving the four
#' peripheral rate constants from one to twelve decades below the
#' absorption rate, which is what a fit does to a compartment the data
#' do not support, leaves the worst error at 6.9e-11 and does not
#' degrade at all.
#'
#' **Beyond an eight-decade spread that includes the FAST rates, both
#' the value and the gradient go wrong.** There the trajectory can be
#' wrong by up to 1.3e-04 of its own maximum against a 300-bit
#' reference, with nothing warning, and the gradient can be `NaN`. That
#' region needs a terminal half-life of seconds beside an absorption
#' rate of hundreds per hour, which is not a drug; and [frm_ode()] is
#' not a usable reference there either, because DLSODA reaches its step
#' limit and disagrees with itself by 1e+05 and more.
#'
#' @section Boundary:
#' A three-compartment model returns `NaN` in the GRADIENT, though not
#' in the value, when the cubic's two smaller roots have lost their
#' precision against each other. `lincmt_disp()` solves the cubic
#' through `acos()`, and a double root puts its argument at exactly
#' one, where `acos()` has an infinite derivative.
#'
#' **The derivative itself is finite there.** What `frm_lincmt()`
#' returns is a symmetric function of the three roots, hence a function
#' of the characteristic polynomial's coefficients, which are
#' polynomials in the rate constants, so the residue formula's poles at
#' a double root are removable and the trajectory is analytic in the
#' rate constants AT the tangency. Measured at an exact double root
#' (`k13` at 1e-300, `k31` on the slow root of the reduced quadratic):
#' `d/d(log ke)` is -9.2779361, stable to eight significant digits at
#' step sizes 1e-4 through 1e-7 and equal to what the tape itself
#' returns one part in 1e8 away. The `NaN` is the implementation
#' routing through the eigenvalues, not a derivative that does not
#' exist.
#'
#' **How often it happens.** Not a knife edge. Over 4000 random draws
#' per box: 0 on a pharmacokinetic box of 1e-2 to 10 per hour, 0 on
#' 1e-4 to 20, 38 on 1e-8 to 50, and 749 on 1e-20 to 50. By spread,
#' 3000 draws per cell: 0 at six decades or fewer, 0.07 percent at
#' eight decades and 1.33 percent at ten, which is one draw in 75.
#'
#' **And the `NaN` is the loud end of a soft region rather than a
#' boundary guarding a right answer.** On ten-decade draws where the
#' gradient comes back FINITE, so nothing fires at all, the value is
#' wrong by up to 1.3e-04 of the trajectory's own maximum. A silent
#' value error outranks a loud `NaN`, and it is the reason this is
#' stated as a domain above rather than patched.
#'
#' **Why it is not clamped.** Holding the `acos()` argument below
#' `1 - 1e-14` makes the gradient finite and right at a genuine double
#' root and costs nothing there (value error 5.177e-13 against the
#' shipped 5.181e-13, and no change at all over 400 ordinary draws).
#' But on 40 draws from the wide region it changes the VALUE by a
#' factor of 2.8e+02 to 2.5e+06. It is the right fix for a real
#' collision and destructive where the roots were already garbage, and
#' telling the two apart needs a test on whether the cubic's
#' coefficients still carry relative precision, which is a comparison
#' on an automatic-differentiation value and is the constraint this
#' whole function is written around. So it is left, and said.
#'
#' **What a user sees.** An optimizer cannot walk INTO the region:
#' fitting three compartments to two-compartment data from six starting
#' points, including `log(k13)` at -20, -30, -35 and -700, every fit
#' returned normally and `nlminb` left `log(k13)` exactly where it
#' began, because the gradient in that direction is zero once `k13`
#' stops mattering. What a user can do is START inside it, and then the
#' fit stops with `nlminb`'s `NA/NaN gradient evaluation`, which is
#' loud and names neither this function nor the cause.
#'
#' A three-compartment model whose five disposition rate constants have
#' all underflowed to exactly zero returns `NaN`, in the value and in
#' the gradient, because the cubic is then a triple root at zero and
#' its trigonometric solution reads `0 / 0`. It fails loudly rather
#' than quietly, and reaching it needs five log rates below about -745
#' at once. It is not offset, because offsetting it by this file's own
#' 1e-150 returns a finite answer at eight-ninths of the true height, a
#' silent factor of 1.125.
#'
#' @param parms A named list of parameters. See "Parameters".
#'   Inside a `bf(nl = TRUE)` body a bare NAME is a request for a
#'   column of the data. The arguments that are not columns
#'   (`ncmt`, `depot`, `output`, `n_ss`) therefore have to be written
#'   as literals there rather than held in a variable.
#' @param times Observation times, one per row of the data.
#' @param group Grouping vector naming the unit that owns one system,
#'   one value per row. `NULL` treats the whole data as one group.
#' @param ncmt Number of compartments in the disposition model: 1, 2 or
#'   3. The depot, if any, is not counted.
#' @param depot Whether a depot compartment feeds the central one by
#'   first-order absorption at rate `ka`.
#' @param output What to return: `"conc"` (the default) is the central
#'   amount divided by `V`, `"central"` and `"depot"` are amounts. One
#'   column only; [frm_ode()]'s per-row `output` is refused by name.
#' @param t0 Initial time, a scalar or one value per row (constant
#'   within group). Every observation time must be at or after it.
#' @param init Amounts present at `t0`, as a named list with elements
#'   `depot` and `central`, each one value per observation (constant
#'   within group) or a single value. `NULL` starts empty.
#' @param events Optional dosing table, exactly [frm_ode()]'s. See
#'   "Dosing".
#' @param event_scale A multiplier on every `events$value`, one value
#'   per observation (constant within group) or a single shared value.
#' @param n_ss How many dosing cycles a steady-state record carries.
#'   `Inf`, the default, is the exact geometric limit and costs one
#'   term. A whole number writes out that many cycles, which is what
#'   [frm_ode()] simulates.
#' @param tv,tv_break Refused, and present so that the refusal names
#'   them. See "What is refused, and why".
#'
#' @return A numeric vector of length `nrow(data)`. On the automatic
#'   differentiation tape it carries the `advector` class.
#'
#' @seealso [frm_ode()] for a system with no closed form, and for every
#'   feature named in "What is refused, and why".
#'
#' @examples
#' # One-compartment oral pharmacokinetics, 100 into the depot every 12
#' # hours, and the same schedule already at steady state.
#' doses <- data.frame(time = 0, state = "depot", value = 100,
#'                     ii = 12, addl = 3L)
#' frm_lincmt(parms = list(ka = 1, ke = 0.2, V = 10),
#'            times = c(6, 18, 30, 42), ncmt = 1, depot = TRUE,
#'            events = doses)
#'
#' frm_lincmt(parms = list(ka = 1, ke = 0.2, V = 10),
#'            times = c(6, 18, 30, 42), ncmt = 1, depot = TRUE,
#'            events = data.frame(time = 0, state = "depot", value = 100,
#'                                ii = 12, ss = TRUE))
#'
#' # The flip-flop case, where the textbook form returns NaN
#' frm_lincmt(parms = list(ka = 0.2, ke = 0.2, V = 10), times = 1:4,
#'            ncmt = 1, depot = TRUE, init = list(depot = 100))
#'
#' # In a formula, with between-subject variability on both rates
#' set.seed(2026)
#' tt <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
#' dd <- data.frame(id = factor(rep(1:6, each = length(tt))),
#'                  time = rep(tt, 6), dose = 100)
#' ka <- exp(rnorm(6, 0, 0.3))[as.integer(dd$id)]
#' ke <- exp(rnorm(6, log(0.2), 0.25))[as.integer(dd$id)]
#' dd$conc <- 100 * ka / (10 * (ka - ke)) *
#'   (exp(-ke * dd$time) - exp(-ka * dd$time)) + rnorm(nrow(dd), 0, 0.3)
#' \donttest{
#' fit <- frm(
#'   bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
#'                                     V = exp(lV)),
#'                        times = time, group = id, ncmt = 1,
#'                        depot = TRUE, init = list(depot = dose)),
#'      lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
#'     gaussian(),
#'   data = dd, start = list(beta = c(0, log(0.25), log(8))))
#' fixef(fit)
#' }
#' @export
frm_lincmt <- function(parms, times, group = NULL, ncmt = 1L,
                       depot = FALSE, output = c("conc", "central",
                                                 "depot"),
                       t0 = 0, init = NULL, events = NULL,
                       event_scale = 1, n_ss = Inf,
                       tv = NULL, tv_break = NULL) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")

  if (!is.null(tv) || !is.null(tv_break)) {
    stop("`tv` is not available in frm_lincmt(). A dynamics input that ",
         "changes with time makes the system time-varying, and the ",
         "closed form is a superposition over the WHOLE history of the ",
         "group, so a rate constant that changed part way through it ",
         "would be applied to doses given before the change. Use ",
         "frm_ode(), which splits the solve at each change point",
         call. = FALSE)
  }
  if (length(ncmt) != 1L || !is.numeric(ncmt) || is.na(ncmt) ||
        !ncmt %in% 1:3) {
    stop("`ncmt` must be 1, 2 or 3: it is how many compartments the ",
         "disposition model has, not counting the depot", call. = FALSE)
  }
  ncmt <- as.integer(ncmt)
  if (length(depot) != 1L || !is.logical(depot) || is.na(depot)) {
    stop("`depot` must be TRUE or FALSE", call. = FALSE)
  }
  if (length(output) > 1L && !identical(output, eval(formals()$output))) {
    stop("`output` names ", length(output), " things. frm_lincmt() ",
         "returns ONE column: \"conc\", \"central\" or \"depot\". ",
         "frm_ode()'s per-row `output`, which reads a different state ",
         "on each row over one shared solve, is not available here; ",
         "call frm_lincmt() once per compartment and combine the ",
         "results in the nonlinear body", call. = FALSE)
  }
  output <- match.arg(output)
  if (identical(output, "depot") && !depot) {
    stop("`output` is \"depot\" but `depot` is FALSE, so there is no ",
         "depot compartment to read", call. = FALSE)
  }
  if (length(n_ss) != 1L || !is.numeric(n_ss) || is.na(n_ss) ||
        n_ss < 1 || (is.finite(n_ss) && n_ss != trunc(n_ss))) {
    stop("`n_ss` must be Inf, or one whole number at least 1: it is ",
         "how many dosing cycles a steady-state record carries",
         call. = FALSE)
  }
  for (nm in c("times", "group", "t0")) {
    if (inherits(get(nm), "advector")) {
      stop("`", nm, "` is an estimated quantity. frm_lincmt() needs it ",
           "as data: it fixes which doses precede which observation, ",
           "which is settled before the tape is built", call. = FALSE)
    }
  }

  times <- as.numeric(times)
  n_obs <- length(times)
  if (!n_obs) stop("`times` is empty", call. = FALSE)
  if (anyNA(times)) stop("`times` contains NA", call. = FALSE)

  if (is.null(group)) {
    gi <- rep(1L, n_obs)
    glab <- "1"
  } else {
    if (length(group) != n_obs) {
      stop("`group` has length ", length(group), " but `times` has ",
           n_obs, call. = FALSE)
    }
    gf <- if (is.factor(group)) droplevels(group) else factor(group)
    gi <- as.integer(gf)
    glab <- levels(gf)
  }
  groups <- split(seq_len(n_obs), gi)
  labels <- glab[as.integer(names(groups))]
  n_grp <- length(groups)
  # the group's first row BY TIME, which is the row frm_ode() reads its
  # dynamics inputs off; they have to be constant within the group
  # anyway, and matching the rule exactly costs one sort
  first <- vapply(groups, function(i) i[order(times[i])][[1L]], 1L)

  n_dep <- as.integer(depot)
  n_state <- ncmt + n_dep
  states <- c(if (depot) "depot", "central",
              if (ncmt >= 2L) "peripheral1", if (ncmt == 3L)
                "peripheral2")

  want_v <- identical(output, "conc")
  pl <- lincmt_rates(parms, ncmt, depot, want_v)
  pcols <- lapply(names(pl), function(nm)
    lincmt_col(pl[[nm]], n_obs, "parms", nm))
  names(pcols) <- names(pl)
  ode_check_constant(pcols, groups, "parms", labels, "frm_lincmt()")
  t0_cols <- ode_columns(t0, n_obs, "t0")
  if (length(t0_cols) != 1L) {
    stop("`t0` must be a single column", call. = FALSE)
  }
  ode_check_constant(t0_cols, groups, "t0", labels, "frm_lincmt()")

  # one value per group for every parameter, in one vector, so that the
  # whole model is evaluated in one vectorized pass rather than in a
  # loop over groups
  per_group <- function(x) if (length(x) == 1L) rep(x, n_grp) else
    x[first]
  pg <- lapply(pcols, per_group)
  rate <- if ("CL" %in% names(pg)) {
    r <- list(ke = pg[["CL"]] / pg[["V"]])
    if (ncmt >= 2L) {
      r[["k12"]] <- pg[["Q2"]] / pg[["V"]]
      r[["k21"]] <- pg[["Q2"]] / pg[["V2"]]
    }
    if (ncmt == 3L) {
      r[["k13"]] <- pg[["Q3"]] / pg[["V"]]
      r[["k31"]] <- pg[["Q3"]] / pg[["V3"]]
    }
    if (depot) r[["ka"]] <- pg[["ka"]]
    r[["V"]] <- pg[["V"]]
    r
  } else pg
  zero_g <- rep(0, n_grp)
  disp <- lincmt_disp(ncmt, rate[["ke"]],
                      if (ncmt >= 2L) rate[["k12"]] else zero_g,
                      if (ncmt >= 2L) rate[["k21"]] else zero_g,
                      if (ncmt == 3L) rate[["k13"]] else zero_g,
                      if (ncmt == 3L) rate[["k31"]] else zero_g)
  ka_g <- if (depot) rate[["ka"]] else zero_g

  init_g <- NULL
  if (!is.null(init)) {
    if (!is.list(init) || inherits(init, "advector") ||
          is.null(names(init)) || any(!nzchar(names(init)))) {
      stop("`init` must be a named list, for example ",
           "init = list(depot = dose). The names are ",
           paste(c(if (depot) "depot", "central"), collapse = " and "),
           call. = FALSE)
    }
    bad <- setdiff(names(init), c(if (depot) "depot", "central"))
    if (length(bad)) {
      stop("`init` names ", paste(bad, collapse = ", "),
           ", and this model's are ",
           paste(c(if (depot) "depot", "central"), collapse = " and "),
           ". frm_lincmt() starts from the ",
           if (depot) "depot or the central compartment" else
             "central compartment",
           " only: an amount in a PERIPHERAL compartment at t0 needs ",
           "that compartment's impulse response, which is the one piece ",
           "of the closed form that has no cancellation-free spelling, ",
           "and a depot needs depot = TRUE. Use frm_ode()",
           call. = FALSE)
    }
    icols <- lapply(names(init), function(nm)
      lincmt_col(init[[nm]], n_obs, "init", nm))
    names(icols) <- names(init)
    ode_check_constant(icols, groups, "init", labels, "frm_lincmt()")
    init_g <- lapply(icols, per_group)
  }

  scale_g <- NULL
  if (!identical(event_scale, 1)) {
    if (is.null(events)) {
      stop("`event_scale` was given but `events` was not; there is ",
           "nothing to scale", call. = FALSE)
    }
    sc <- ode_columns(event_scale, n_obs, "event_scale")
    if (length(sc) != 1L) {
      stop("`event_scale` must be a single column", call. = FALSE)
    }
    ode_check_constant(sc, groups, "event_scale", labels,
                       "frm_lincmt()")
    scale_g <- per_group(sc[[1L]])
  }

  ev_by_group <- if (is.null(events)) NULL else
    ode_split_events(events, labels, n_state, states)

  # ------------------------------------------------------------------
  # The term table. Everything in it is data: which doses precede which
  # observation, and by how long, is settled here, once.
  per_row <- vector("list", n_obs)
  for (g in seq_len(n_grp)) {
    idx <- groups[[g]]
    tstart <- as.numeric(if (length(t0_cols[[1L]]) == 1L) t0_cols[[1L]]
                         else t0_cols[[1L]][first[[g]]])
    if (min(times[idx]) < tstart) {
      stop("group '", labels[[g]], "' has an observation time (",
           format(min(times[idx])), ") before t0 (", format(tstart),
           "); frm_lincmt() runs forward from t0", call. = FALSE)
    }
    ev <- if (is.null(ev_by_group)) NULL else ev_by_group[[labels[[g]]]]
    if (!is.null(ev) && !nrow(ev)) ev <- NULL
    if (!is.null(ev)) {
      if (any(ev[["time"]] < tstart)) {
        stop("group '", labels[[g]], "' has an event at time ",
             format(min(ev[["time"]])), ", before t0 (", format(tstart),
             "); frm_lincmt() runs forward from t0", call. = FALSE)
      }
      lincmt_check_events(ev, labels[[g]], n_dep, n_state)
    }
    rst <- if (is.null(ev)) numeric(0) else
      ev[["time"]][ev[["method"]] == "reset"]
    sst <- if (is.null(ev)) numeric(0) else ev[["time"]][ev[["ss"]]]
    dose <- if (is.null(ev)) NULL else
      ev[ev[["method"]] == "add", , drop = FALSE]
    # an infusion that spans a restart would have to be carried across
    # it, and superposition drops everything before a restart
    for (b in c(rst, sst)) {
      if (!is.null(dose) && any(dose[["duration"]] > 0 &
                                dose[["time"]] < b &
                                dose[["time"]] + dose[["duration"]] > b)) {
        stop("group '", labels[[g]], "' has an infusion running at ",
             "time ", format(b), ", where the schedule restarts (a ",
             "\"reset\" row, or a steady-state row). frm_lincmt() drops ",
             "every dose before a restart, so an infusion that spans ",
             "one would lose the part still to be delivered. Split the ",
             "infusion at the restart, or use frm_ode()",
             call. = FALSE)
      }
    }
    # 1 is the depot and 2 the central compartment, whatever position
    # the schedule numbered them at
    cmt_of <- function(s) if (depot) ifelse(s == 1L, 1L, 2L) else
      rep(2L, length(s))
    for (i in idx) {
      tt <- times[[i]]
      # An observation at a steady-state time reads the run-in's own
      # trough, and one at a reset time reads the state BEFORE the
      # reset: that asymmetry is frm_ode()'s, and it is matched rather
      # than tidied. A reset at the same instant as a steady-state row
      # is applied after it, so the reset wins.
      r_ss <- if (length(sst) && any(sst <= tt)) max(sst[sst <= tt]) else
        -Inf
      r_rs <- if (length(rst) && any(rst < tt)) max(rst[rst < tt]) else
        -Inf
      r <- max(tstart, r_ss, r_rs)
      is_ss <- r_ss == r && r_rs < r_ss
      use_init <- r == tstart && !is_ss && !isTRUE(r_rs == tstart) &&
        !is.null(init_g)
      keep <- if (is.null(dose)) integer(0) else
        which(dose[["time"]] >= r & dose[["time"]] < tt)
      u <- numeric(0); cm <- integer(0); du <- numeric(0)
      ii <- numeric(0); kd <- integer(0); am <- numeric(0)
      sl <- logical(0)
      if (length(keep)) {
        u <- tt - dose[["time"]][keep]
        cm <- cmt_of(dose[["state"]][keep])
        du <- dose[["duration"]][keep]
        ii <- rep(0, length(keep))
        kd <- rep(1L, length(keep))
        am <- dose[["value"]][keep]
        sl <- rep(TRUE, length(keep))
      }
      if (is_ss) {
        h <- which(ev[["ss"]] & ev[["time"]] == r)[[1L]]
        sii <- ev[["ii"]][[h]]
        scm <- cmt_of(ev[["state"]][[h]])
        if (is.finite(n_ss)) {
          k <- seq_len(as.integer(n_ss))
          u <- c(u, tt - r + k * sii)
          cm <- c(cm, rep(scm, length(k)))
          du <- c(du, rep(ev[["duration"]][[h]], length(k)))
          ii <- c(ii, rep(0, length(k)))
          kd <- c(kd, rep(1L, length(k)))
          am <- c(am, rep(ev[["value"]][[h]], length(k)))
          sl <- c(sl, rep(TRUE, length(k)))
        } else {
          u <- c(u, tt - r + sii)
          cm <- c(cm, scm)
          du <- c(du, ev[["duration"]][[h]])
          ii <- c(ii, sii)
          kd <- c(kd, 2L)
          am <- c(am, ev[["value"]][[h]])
          sl <- c(sl, TRUE)
        }
      } else if (use_init) {
        for (nm in names(init_g)) {
          u <- c(u, tt - tstart)
          cm <- c(cm, if (identical(nm, "depot")) 1L else 2L)
          du <- c(du, 0)
          ii <- c(ii, 0)
          kd <- c(kd, 1L)
          am <- c(am, NA_real_)
          sl <- c(sl, FALSE)
        }
      }
      per_row[[i]] <- list(u = u, cmt = cm, dur = du, ii = ii,
                           kind = kd, amt = am, scaled = sl,
                           init = if (use_init) names(init_g) else
                             character(0))
    }
  }

  nk <- vapply(per_row, function(z) length(z[["u"]]), 1L)
  kmax <- max(1L, max(nk))
  n_cell <- n_obs * kmax
  # padding is a dose of zero into the central compartment at lag zero:
  # it evaluates without a branch and contributes nothing
  cu <- numeric(n_cell)
  ccm <- rep(2L, n_cell)
  cdu <- numeric(n_cell)
  cii <- numeric(n_cell)
  ckd <- rep(1L, n_cell)
  camt <- numeric(n_cell)
  cscl <- rep(FALSE, n_cell)
  cini <- rep(NA_character_, n_cell)
  cg <- rep(gi, kmax)
  for (i in seq_len(n_obs)) {
    z <- per_row[[i]]
    if (!nk[[i]]) next
    at <- i + (seq_len(nk[[i]]) - 1L) * n_obs
    cu[at] <- z[["u"]]
    ccm[at] <- z[["cmt"]]
    cdu[at] <- z[["dur"]]
    cii[at] <- z[["ii"]]
    ckd[at] <- z[["kind"]]
    camt[at] <- ifelse(is.na(z[["amt"]]), 0, z[["amt"]])
    cscl[at] <- z[["scaled"]]
    if (length(z[["init"]])) {
      cini[at[is.na(z[["amt"]])]] <- z[["init"]]
    }
  }

  # ------------------------------------------------------------------
  # One vectorized pass.
  out_idx <- if (identical(output, "depot")) 1L else 2L
  resp <- lincmt_resp(cu, cdu, cii, ckd, ccm,
                      lapply(disp[["lam"]], function(x) x[cg]),
                      lapply(disp[["coef"]], function(x)
                        if (length(x) == 1L) x else x[cg]),
                      if (depot) ka_g[cg] else 0, out_idx)
  amt <- if (is.null(scale_g)) camt else {
    w <- scale_g[cg]
    w[!cscl] <- 1
    camt * w
  }
  if (!is.null(init_g)) {
    for (nm in names(init_g)) {
      at <- which(cini == nm)
      if (length(at)) amt[at] <- init_g[[nm]][cg[at]]
    }
  }
  cell <- amt * resp
  res <- if (kmax == 1L) cell else
    RTMB::rowSums(RTMB::matrix(cell, n_obs, kmax))
  if (identical(output, "conc")) res <- res / rate[["V"]][gi]
  res
}
