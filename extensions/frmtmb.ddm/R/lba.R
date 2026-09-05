#' The linear ballistic accumulator
#'
#' A race between `n` accumulators, for choices with more than two
#' alternatives. Each accumulator rises in a straight line, with no
#' within-trial noise, from a start point drawn uniformly on `(0, A)`
#' toward a common threshold `b`. Its rate is drawn once per trial from
#' a normal distribution. The first accumulator to reach the threshold
#' is the response, and the observed time is its arrival time plus a
#' non-decision time.
#'
#' Because a trial's outcome is decided by one draw per accumulator
#' rather than by a path, the likelihood is closed form for any `n`.
#' That is what the family is for: [wiener()] is a diffusion between two
#' absorbing boundaries and so admits exactly two responses, while this
#' one takes as many as the design has.
#'
#' @section The model:
#' Write the decision time as \eqn{t = rt - ndt}. For one accumulator
#' with drift mean \eqn{v} and drift standard deviation \eqn{s}, put
#'
#' \deqn{g = \frac{b - v t}{s t}, \qquad
#'       h = \frac{b - A - v t}{s t},}
#'
#' so that \eqn{g} and \eqn{h} are the standardized drifts that place the
#' accumulator exactly at the threshold and exactly at the top of the
#' start-point range at time \eqn{t}. The defective density and
#' distribution function of its arrival time are
#'
#' \deqn{f(t) = \frac{v\,(\Phi(g) - \Phi(h)) +
#'                    s\,(\phi(h) - \phi(g))}{A\,c},}
#' \deqn{S(t) = \frac{\Phi(g) - \Phi(-v/s) +
#'                    \frac{s t}{A}\,(h\,(\Phi(g) - \Phi(h)) -
#'                    (\phi(h) - \phi(g)))}{c},}
#'
#' where \eqn{S} is the probability of arriving after \eqn{t} and
#' \eqn{c} is the drift normalizing constant of the next section. A
#' trial on which accumulator \eqn{j} responded at time \eqn{rt}
#' contributes
#'
#' \deqn{\log f_j(t) + \sum_{i \neq j} \log S_i(t),}
#'
#' the density that the observed accumulator arrived then, times the
#' probability that none of the others had arrived yet.
#'
#' @section Identification:
#' The model is invariant to a common rescaling: multiplying `A`, `b`,
#' every drift mean and the drift standard deviation by one constant
#' leaves the distribution of `(choice, rt)` unchanged. One quantity
#' must therefore be fixed, and this family fixes the drift standard
#' deviation, which is the usual choice.
#'
#' `sd_v` is a family argument rather than a distributional parameter,
#' so the fixed quantity is written where the model is written and is
#' carried on the family object, not buried in a default. Change it and
#' every other parameter moves in proportion; it is a choice of units.
#'
#' The other convention in the literature constrains the drift means to
#' sum to one. This family does not offer it, because it would make each
#' accumulator's drift a function of every other accumulator's linear
#' predictor, and the whole point of the family is that each drift takes
#' its own formula.
#'
#' @section Drift rates below zero:
#' An accumulator whose drift is not positive never reaches the
#' threshold. Implementations differ on what to do about that. This
#' family follows `rtdists`, the canonical R implementation, and its
#' default `posdrift = TRUE`: the drift distribution is a normal
#' truncated to positive values, so \eqn{c = \Phi(v/s)} divides both
#' \eqn{f} and \eqn{S}. Every accumulator then arrives eventually, the
#' choice probabilities sum to one, and the likelihood is proper.
#'
#' `posdrift = FALSE` leaves the drift distribution untruncated, so
#' \eqn{c = 1} and an accumulator has probability \eqn{\Phi(-v/s)} of
#' never arriving. The response distribution is then defective: the
#' choice probabilities sum to less than one, the missing mass being
#' trials on which no accumulator would ever respond.
#'
#' The two conventions differ by more than a constant, because the
#' constant is a different one for each accumulator and depends on that
#' accumulator's own drift. A fit under one convention is not a
#' reparameterization of a fit under the other. If you are comparing
#' against another package, check which it uses before comparing
#' numbers: `rtdists` truncates by default, and so does this family.
#'
#' @section Parameters:
#' \describe{
#'   \item{`v1`, ..., `vn`}{Drift means, one per accumulator. Identity
#'     link, because a drift mean is signed: it is the mean of the
#'     normal distribution before truncation, so a negative value is
#'     meaningful and describes an accumulator that rarely wins. These
#'     are the primary parameters, so the main formula goes to all of
#'     them and each gets its own coefficients. Give one its own
#'     formula to move it alone.}
#'   \item{`A`}{Upper end of the start-point range, so start points are
#'     uniform on `(0, A)`. Log link.}
#'   \item{`k`}{Distance from the top of the start-point range to the
#'     threshold, so that `b = A + k`. Log link.}
#'   \item{`ndt`}{Non-decision time. Bounded link; see below.}
#' }
#'
#' The threshold is parameterized as `A + k` rather than directly,
#' which is the spelling the LBA literature uses and which makes
#' `b > A` structural: a log link keeps `k` positive at every value of
#' its linear predictor, so a threshold inside the start-point range,
#' where a fraction of trials would start already finished, is not a
#' state the optimizer can reach. Pinning `k` to zero or a negative
#' constant is refused when the constant is checked against the log
#' link's range.
#'
#' @section Non-decision time:
#' The density is zero at and below `ndt`, so the likelihood has a hard
#' edge at `ndt = min(rt)`. As [wiener()] does, `ndt` gets a logit
#' scaled onto `(0, max_ndt)` instead of a log link, which makes the
#' constraint structural. `max_ndt` defaults to the smallest observed
#' response time, taken when the model frame is assembled. Pass it
#' explicitly to pin the bound, which matters if you will `predict()`
#' on new data whose minimum differs.
#'
#' @section The response, and why not `dec()`:
#' A trial is a `(choice, time)` pair. The time is the response and the
#' choice is per-row data, which reaches the family through `vint()`:
#'
#' ```
#' frm(bf(rt | vint(choice) ~ cond), family = lba(3), data = dat)
#' ```
#'
#' `vint1` is the ACCUMULATOR INDEX for this family: a whole number in
#' `1..n` naming which accumulator reached the threshold, counting from
#' one. A factor is not accepted, so recode it with
#' `as.integer(factor(choice))` and check that the level order matches
#' the accumulator numbering. Omitting `vint()` is refused by name,
#' because the density indexes it.
#'
#' `dec()` does not apply here. That addition term is the two-boundary
#' families' spelling: it carries a 0/1 indicator saying which of two
#' boundaries a trial ended at, and a race has no boundaries and no
#' fixed count of two. Passing `dec()` to this family is refused,
#' because a decision indicator has two levels and a race of `n`
#' accumulators needs `n`. The two are not different spellings of one
#' idea: `dec()` names a side, `vint1` here names a winner. Read a
#' `dec()` example from [wiener()] as 0/1, and this one as `1..n`.
#'
#' @section Accuracy:
#' The pieces the likelihood is built from are each written in the form
#' that keeps its digits, which is not always the form the papers print.
#'
#' The single-accumulator density agrees with `rtdists::dlba_norm()` to
#' better than 1e-11 relative wherever `rtdists` itself is accurate,
#' which is decision times long enough to keep `h` below about
#' 4.5. Past that the two DIVERGE BY DESIGN and this one is the
#' better of them, adjudicated at 200 bits: already at `h` = 7.6,
#' inside the band the literature would call safe, `rtdists` is
#' 1.3e-4 wrong where this family is 4.6e-15 wrong.
#' `rtdists` writes the normal difference `Phi(g) - Phi(h)` as a
#' subtraction of two lower tails, which returns exactly zero once
#' `Phi(h)` rounds to one; this family writes it in log space, where
#' nothing saturates. Against a 200-bit reference on a grid of 1144
#' points, the subtractive form is wrong by a factor of one on 68 of
#' the 80 rows with `h` at or above 8, returning exactly zero, and is
#' already 2.7e-3 wrong on 8 rows well inside the supposedly safe
#' region. The log-space form holds 5.0e-14 in the bulk and 3.6e-4 in
#' the tail, with no exact zeros.
#'
#' Four of those 1144 points still exceed 1e-6 relative, the worst at
#' 3.6e-4, so this is a large improvement rather than a proof.
#'
#' The survival function is likewise not `1 - rtdists::plba_norm()`:
#' that difference loses all its digits once a fast competitor has
#' almost certainly finished, returning exactly zero where the survival
#' is around 1e-19, which would send a loser's log-contribution to
#' `-Inf` on ordinary data. It is written out directly instead, and
#' agrees with a high-accuracy quadrature of the density to better than
#' 1e-9 relative down to survivals of 1e-23.
#'
#' Why this was worth doing rather than documenting as a limit: the
#' subtractive form's error is an error in the VALUE alone. A tape
#' differentiates the function that was written, not the one that was
#' computed, so where the value saturates the gradient stops describing
#' the surface the value traces. That mismatch does not move the point
#' estimate, which the likelihood keeps out of the affected region by
#' collapsing the fastest row's density, but it does reach the Hessian,
#' and standard errors were up to 16.5 percent off toward the
#' non-decision-time bound before the change.
#'
#' Below all of this the density underflows to zero in double
#' precision; the log-density is floored rather than returning `-Inf`,
#' so the optimizer sees a finite wall instead of a hole. Decision
#' times at or below zero, which `predict()` on faster new data can
#' reach, are floored the same way.
#'
#' @param n Number of accumulators, so the number of response
#'   alternatives. At least 2.
#' @param sd_v Fixed drift standard deviation, the quantity that
#'   identifies the scale. One positive number, or one per accumulator.
#' @param posdrift Truncate the drift distribution at zero? `TRUE`, the
#'   default, matches `rtdists` and makes the likelihood proper.
#' @param max_ndt Upper bound for the non-decision time, in the units of
#'   the response. `NULL`, the default, takes it from the data.
#'
#' @return A `frmtmb_family`.
#'
#' @references
#' Brown, S. D. and Heathcote, A. (2008). The simplest complete model of
#' choice response time: Linear ballistic accumulation. *Cognitive
#' Psychology*, 57(3), 153-178.
#'
#' Donkin, C., Brown, S. D. and Heathcote, A. (2009). The overconstraint
#' of response time models: Rethinking the scaling problem.
#' *Psychonomic Bulletin & Review*, 16(6), 1129-1135.
#'
#' @seealso [lba_simulate()] to generate from the model, and [wiener()]
#'   for the two-choice diffusion.
#'
#' @examples
#' set.seed(1)
#' dat <- lba_simulate(400, v = c(2.4, 1.6, 1.0), A = 0.5, k = 0.4,
#'                     ndt = 0.2)
#' fit <- frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = dat)
#' fixef(fit)
#'
#' @export
lba <- function(n, sd_v = 1, posdrift = TRUE, max_ndt = NULL) {
  if (missing(n) || !is.numeric(n) || length(n) != 1L || is.na(n) ||
      n != round(n) || n < 2) {
    stop("lba(): `n` is the number of accumulators, one whole number ",
         "of 2 or more, e.g. lba(3) for a three-alternative choice.",
         call. = FALSE)
  }
  n <- as.integer(n)
  if (!is.numeric(sd_v) || !all(is.finite(sd_v)) || any(sd_v <= 0) ||
      !length(sd_v) %in% c(1L, n)) {
    stop("lba(): `sd_v` fixes the scale of the model and must be one ",
         "positive finite number, or one per accumulator (", n, " of ",
         "them). It is not estimated: something has to be held fixed ",
         "or the model is not identified.", call. = FALSE)
  }
  sd_v <- rep(as.numeric(sd_v), length.out = n)
  if (!is.logical(posdrift) || length(posdrift) != 1L || is.na(posdrift)) {
    stop("lba(): `posdrift` says whether drift rates are truncated at ",
         "zero and must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(max_ndt)) {
    if (!is.numeric(max_ndt) || length(max_ndt) != 1L ||
        !is.finite(max_ndt) || max_ndt <= 0) {
      stop("lba(): `max_ndt` bounds the non-decision time and must be ",
           "one positive finite number, or NULL to read it off the ",
           "response.", call. = FALSE)
    }
  }

  vp <- paste0("v", seq_len(n))
  dpn <- c(vp, "A", "k", "ndt")
  lk <- c(rep(list("identity"), n), list("log"), list("log"), list("log"))
  names(lk) <- dpn

  fam <- frmtmb::custom_family(
    "lba",
    dpars = dpn,
    links = lk,
    lpdf = function(y, dpars, aterms) {
      lba_race_lpdf(y - dpars[["ndt"]], aterms[["vint1"]],
                    lba_law, lba_pars(dpars, vp, sd_v, posdrift))
    },
    valid_y = function(y, aterms) lba_check_response(y, aterms, n),
    family_finalize = function(fam, y, aterms) {
      lba_finalize(fam, y, max_ndt)
    },
    required_aterms = "vint1",
    init_dpars = lba_inits(vp),
    type = "continuous",
    sim = function(dpars, aterms, n_) {
      lba_sim_rt(dpars, aterms, n_, vp, sd_v, posdrift)
    },
    primary_dpars = vp)

  # Carried so that a reader of the fitted object can see what was held
  # fixed, and so that the simulator and the post-fit helpers do not
  # have to be told again.
  fam[["lba_n"]] <- n
  fam[["lba_sd_v"]] <- sd_v
  fam[["lba_posdrift"]] <- posdrift
  fam
}

#' The one accumulator law the race is written against.
#'
#' `ldens` and `lsurv` are all `lba_race_lpdf()` knows about, so a
#' different single-accumulator law with the same two functions races
#' without touching the race itself. Each takes a decision time and one
#' accumulator's parameter list and returns a vector the same length.
#'
#' @noRd
lba_law <- list(
  ldens = function(t, p) {
    # RTMB's dnorm, not stats', because this runs on the tape: the stats
    # version takes the AD class off and the derivative with it.
    phi <- RTMB::dnorm
    g <- (p$b - p$v * t) / (p$s * t)
    h <- (p$b - p$A - p$v * t) / (p$s * t)
    dPhi <- lba_phidiff(g, h)
    dphi <- phi(h) - phi(g)
    f <- (p$v * dPhi + p$s * dphi) / (p$A * lba_denom(p))
    # Floored rather than clamped at zero: the race multiplies this by a
    # 0/1 mask, and 0 * -Inf is NaN, which would take the whole tape with
    # it on the first row whose density underflowed.
    log(lba_atleast(f, 1e-300))
  },
  lsurv = function(t, p) {
    Phi <- RTMB::pnorm
    phi <- RTMB::dnorm
    g <- (p$b - p$v * t) / (p$s * t)
    h <- (p$b - p$A - p$v * t) / (p$s * t)
    dPhi <- lba_phidiff(g, h)
    dphi <- phi(h) - phi(g)
    # Written out rather than as 1 - F. The two agree in the bulk, but
    # F rounds to exactly 1 while the survival is still around 1e-19,
    # and a loser that has probably-but-not-certainly finished is an
    # ordinary row, not an extreme one.
    #
    # `sf` is the mass that arrives strictly after t. An accumulator that
    # never arrives has also not arrived by t, so which of the two ways
    # that mass is accounted for is exactly the drift convention: under
    # truncation there is none of it and sf is renormalized, without it
    # the never-arriving mass Phi(-v/s) is added back.
    q <- Phi(p$v / p$s)
    sf <- lba_phidiff(g, -p$v / p$s) + (p$s * t / p$A) * (h * dPhi - dphi)
    s <- if (p$posdrift) sf / lba_atleast(q, 1e-10) else sf + Phi(-p$v / p$s)
    log(lba_atleast(s, 1e-300))
  })

#' `Phi(hi) - Phi(lo)` for `hi >= lo`, in log space.
#'
#' The published formula, and `rtdists`, spell this as the difference of
#' two LOWER tails. That form loses everything once `lo` passes about
#' 8.3, where `pnorm(lo)` rounds to exactly one: the difference
#' collapses to zero, the density drops its `v (Phi(hi) - Phi(lo))`
#' term, and it is understated by a factor of about `1 + v / (s lo)`.
#' The difference of two UPPER tails is algebraically the same and
#' fails in the mirror image, at `hi` below -8.3.
#'
#' Neither needs to be chosen. Writing the same quantity as
#'
#' ```
#' la <- pnorm(-lo, log.p = TRUE)   # log of the upper tail at lo
#' lb <- pnorm(-hi, log.p = TRUE)   # log of the upper tail at hi
#' exp(la) * -expm1(lb - la)
#' ```
#'
#' keeps full RELATIVE precision at both ends and needs no comparison,
#' no smooth blend and no tape configuration. `pnorm(log.p = TRUE)`
#' does not saturate, because the log of a tail near one is a small
#' negative number a double represents exactly where the tail itself
#' rounds to one; `expm1` then removes the cancellation. The result is
#' a product of two factors each carrying full precision, rather than a
#' difference of two numbers that nearly agree, which is what every
#' form that subtracts eventually founders on.
#'
#' Measured against a 200-bit Rmpfr reference on this package's own
#' density grid (1144 rows with a positive true value): for `lo < 8`
#' the lower-tail form reaches 2.7e-3 relative with 8 rows worse than
#' 1e-6, and this form reaches 5.0e-14 with none; for `lo >= 8` the
#' lower-tail form is wrong by a factor of one on 68 of 80 rows, which
#' return exactly zero, and this form returns no exact zeros and stays
#' within 3.6e-4. The 8.3 threshold understated the old problem: the
#' lower-tail form was already parts-per-thousand wrong well inside it.
#'
#' The caveat, kept because it is real: 4 of those 1144 rows still
#' exceed 1e-6 relative, the worst at 3.6e-4. A large improvement, not
#' a proof of correctness to machine precision.
#'
#' An earlier revision of this file argued that no fix existed, on the
#' grounds that the two subtractive forms fail in mirror-image regimes
#' and RTMB refuses comparison on AD types so neither can be selected.
#' Both halves were true and the conclusion did not follow: it assumed
#' the answer had to come from a subtraction. It does not.
#'
#' @noRd
lba_phidiff <- function(hi, lo) {
  la <- RTMB::pnorm(-lo, log.p = TRUE)
  lb <- RTMB::pnorm(-hi, log.p = TRUE)
  exp(la) * -expm1(lb - la)
}

#' The drift normalizing constant.
#'
#' Under truncation this is the probability that the accumulator's drift
#' is positive, so that it arrives at all. The floor is `rtdists`'s own,
#' kept so that the two agree exactly: below a drift of about -6.4
#' standard deviations it changes the density, and both packages change
#' it the same way.
#'
#' @noRd
lba_denom <- function(p) {
  if (!p$posdrift) return(1)
  lba_atleast(RTMB::pnorm(p$v / p$s), 1e-10)
}

#' `pmax(x, lo)` for a floor `lo` many orders of magnitude below `x`.
#'
#' The obvious spelling, `0.5 * (x + lo + abs(x - lo))`, is wrong here
#' and silently so: with `lo` at 1e-300 and `x` of order one, both
#' `x + lo` and `abs(x - lo)` round to `x` and the floor is annihilated,
#' leaving `pmax(x, 0)`. A negative `x` then floors to zero and its log
#' to `-Inf`, which is the failure the floor existed to prevent. Adding
#' the floor last, outside the cancelling sum, keeps it.
#'
#' @noRd
lba_atleast <- function(x, lo) lo + 0.5 * ((x - lo) + abs(x - lo))

#' Split the distributional parameters into one list per accumulator.
#'
#' @noRd
lba_pars <- function(dpars, vp, sd_v, posdrift) {
  b <- dpars[["A"]] + dpars[["k"]]
  out <- vector("list", length(vp))
  for (j in seq_along(vp)) {
    out[[j]] <- list(v = dpars[[vp[j]]], A = dpars[["A"]], b = b,
                     s = sd_v[j], posdrift = posdrift)
  }
  out
}

#' The race: the winner's log density plus every loser's log survival.
#'
#' The winner is data, so the selection is a numeric 0/1 mask rather
#' than a branch, and the same tape serves every row. Nothing here is
#' specific to the linear ballistic accumulator; `law` is the only
#' place the single-accumulator distribution enters.
#'
#' @noRd
lba_race_lpdf <- function(t, choice, law, accs) {
  # A decision time of zero or less has no accumulator law at all: the
  # standardized drifts change sign and the density is meaningless
  # rather than small. The bounded `ndt` link keeps the fitted model out
  # of that region, but predict() on new data holding the training
  # bound can still reach it, so the time is floored and those rows get
  # the same finite wall an underflowed density gets.
  t <- lba_atleast(t, 1e-12)
  ll <- 0 * t
  for (j in seq_along(accs)) {
    won <- as.numeric(choice == j)
    ll <- ll + won * law$ldens(t, accs[[j]]) +
      (1 - won) * law$lsurv(t, accs[[j]])
  }
  ll
}

#' Response and choice validation.
#'
#' @noRd
lba_check_response <- function(y, aterms, n) {
  if (any(!is.finite(y)) || any(y <= 0)) {
    stop("lba: the response must be a strictly positive, finite ",
         "response time. A time of zero or less has no decision in it ",
         "for any non-decision time.", call. = FALSE)
  }
  ch <- aterms[["vint1"]]
  if (any(!is.finite(ch)) || any(ch != round(ch)) ||
      any(ch < 1) || any(ch > n)) {
    bad <- unique(ch[!is.finite(ch) | ch != round(ch) | ch < 1 | ch > n])
    stop("lba(", n, "): the vint() choice indicator names which ",
         "accumulator responded and must be a whole number from 1 to ",
         n, ". Saw ", paste(utils::head(sort(bad), 5), collapse = ", "),
         ". A factor is not accepted; recode it with ",
         "as.integer(factor(choice)) and check the level order ",
         "matches the accumulator numbering.", call. = FALSE)
  }
  invisible(NULL)
}

#' Fit the non-decision-time link to the response.
#'
#' Uses the `family_finalize` seam rather than writing the bound into an
#' environment the link closures read later: the family object then
#' says what it is, and the bound cannot be read before it is set.
#'
#' @noRd
lba_finalize <- function(fam, y, max_ndt) {
  ub <- if (is.null(max_ndt)) min(y) else max_ndt
  if (!is.null(max_ndt) && ub > min(y)) {
    stop("lba: max_ndt = ", format(ub), " is above the fastest response ",
         "(", format(min(y)), "). No accumulator can arrive before the ",
         "non-decision time, so a bound above the fastest response ",
         "admits parameter values at which that trial has no ",
         "likelihood.", call. = FALSE)
  }
  fam[["links"]][["ndt"]] <- list(
    name = paste0("scaled_logit(0, ", signif(ub, 4), ")"),
    linkfun = function(mu) log(mu / (ub - mu)),
    linkinv = function(eta) ub / (1 + exp(-eta)),
    mu_eta = function(eta) {
      p <- 1 / (1 + exp(-eta))
      ub * p * (1 - p)
    })
  fam
}

#' Starting values.
#'
#' The drifts start apart rather than equal: all-equal drifts are a
#' saddle for the race, since every accumulator then has the same
#' chance and no gradient distinguishes them.
#'
#' @noRd
lba_inits <- function(vp) {
  out <- list(A = function(y, aterms) 0.5,
              k = function(y, aterms) 0.5,
              ndt = function(y, aterms) 0.5 * min(y))
  for (j in seq_along(vp)) {
    out[[vp[j]]] <- local({
      jj <- j
      function(y, aterms) 1.5 - 0.1 * (jj - 1)
    })
  }
  out[c(vp, "A", "k", "ndt")]
}

#' Draw response times conditional on each row's observed choice.
#'
#' The natural draw from this model is a joint `(choice, time)` pair,
#' but `simulate()` replaces the response only and the choice is data,
#' so the draw is conditioned on it by rejection. See [lba_simulate()]
#' for the joint draw, which is what a user usually wants.
#'
#' @noRd
lba_sim_rt <- function(dpars, aterms, n_, vp, sd_v, posdrift) {
  ch <- rep(as.numeric(aterms[["vint1"]]), length.out = n_)
  A <- rep(as.numeric(dpars[["A"]]), length.out = n_)
  b <- A + rep(as.numeric(dpars[["k"]]), length.out = n_)
  ndt <- rep(as.numeric(dpars[["ndt"]]), length.out = n_)
  V <- vapply(vp, function(p) rep(as.numeric(dpars[[p]]), length.out = n_),
              numeric(n_))
  dim(V) <- c(n_, length(vp))
  out <- rep(NA_real_, n_)
  todo <- seq_len(n_)
  for (round in 1:200) {
    if (!length(todo)) break
    dr <- lba_race_draw(A[todo], b[todo], V[todo, , drop = FALSE],
                        sd_v, posdrift)
    hit <- dr$choice == ch[todo]
    out[todo[hit]] <- dr$time[hit] + ndt[todo[hit]]
    todo <- todo[!hit]
  }
  if (length(todo)) {
    stop("lba: could not draw a response time for ", length(todo),
         " row(s) whose observed choice the fitted parameters almost ",
         "never produce. simulate() conditions each draw on that row's ",
         "vint() choice, so a choice the model gives a vanishing ",
         "probability has no draw to give. Use lba_simulate() for an ",
         "unconditional draw of choice and time together.",
         call. = FALSE)
  }
  out
}

#' One unconditional pass of the race, in the generative parameterization.
#'
#' @noRd
lba_race_draw <- function(A, b, V, sd_v, posdrift) {
  m <- nrow(V)
  nacc <- ncol(V)
  tim <- matrix(0, m, nacc)
  for (j in seq_len(nacc)) {
    x <- stats::runif(m, 0, A)
    if (posdrift) {
      lo <- stats::pnorm(0, V[, j], sd_v[j])
      d <- stats::qnorm(lo + stats::runif(m) * (1 - lo), V[, j], sd_v[j])
    } else {
      d <- stats::rnorm(m, V[, j], sd_v[j])
    }
    # A drift at or below zero never arrives; Inf is the honest time and
    # loses every race it is in.
    tim[, j] <- ifelse(d > 0, (b - x) / d, Inf)
  }
  w <- max.col(-tim, ties.method = "first")
  list(choice = w, time = tim[cbind(seq_len(m), w)])
}

#' Simulate from a linear ballistic accumulator
#'
#' Draws from the generative process directly: a uniform start point and
#' a normal drift rate per accumulator, and whichever reaches the
#' threshold first. This is the joint draw of choice and time, which is
#' what the model produces and what a fit needs; `simulate()` on a
#' fitted object instead redraws only the time, holding each row's
#' observed choice.
#'
#' @param n Number of trials.
#' @param v Drift means, one per accumulator. Its length sets the number
#'   of accumulators. May be a matrix with one row per trial, for a
#'   drift that varies with a covariate.
#' @param A Upper end of the start-point range.
#' @param k Distance from the top of the start-point range to the
#'   threshold, so the threshold is `A + k`.
#' @param ndt Non-decision time.
#' @param sd_v Drift standard deviation, fixed. One number or one per
#'   accumulator.
#' @param posdrift Truncate the drift distribution at zero? Matches the
#'   [lba()] argument of the same name.
#'
#' @return A data frame with `choice` (the accumulator that responded,
#'   from 1) and `rt`.
#'
#' @examples
#' set.seed(1)
#' dat <- lba_simulate(500, v = c(2.5, 1.5, 1.0), A = 0.5, k = 0.4,
#'                     ndt = 0.2)
#' table(dat$choice)
#' tapply(dat$rt, dat$choice, mean)
#'
#' @export
lba_simulate <- function(n, v, A = 0.5, k = 0.4, ndt = 0.2, sd_v = 1,
                         posdrift = TRUE) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("lba_simulate(): `n` is the number of trials to draw, one ",
         "whole number of 1 or more.", call. = FALSE)
  }
  n <- as.integer(n)
  V <- if (is.matrix(v)) v else matrix(as.numeric(v), n, length(v),
                                       byrow = TRUE)
  if (nrow(V) != n) {
    stop("lba_simulate(): a matrix `v` gives one row of drift means per ",
         "trial, so it needs ", n, " rows, not ", nrow(V), ".",
         call. = FALSE)
  }
  nacc <- ncol(V)
  if (nacc < 2) {
    stop("lba_simulate(): `v` needs a drift mean for each of at least ",
         "two accumulators; a race of one has nothing to lose to.",
         call. = FALSE)
  }
  if (!length(sd_v) %in% c(1L, nacc) || any(sd_v <= 0)) {
    stop("lba_simulate(): `sd_v` must be one positive number, or one ",
         "for each of the ", nacc, " accumulators.", call. = FALSE)
  }
  sd_v <- rep(as.numeric(sd_v), length.out = nacc)
  A <- rep(as.numeric(A), length.out = n)
  k <- rep(as.numeric(k), length.out = n)
  if (any(A <= 0) || any(k <= 0)) {
    stop("lba_simulate(): `A` and `k` must both be positive, so that ",
         "the threshold A + k lies above the start-point range and no ",
         "trial begins already finished.", call. = FALSE)
  }
  dr <- lba_race_draw(A, A + k, V, sd_v, posdrift)
  data.frame(choice = dr$choice,
             rt = dr$time + rep(as.numeric(ndt), length.out = n))
}
