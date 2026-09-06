#' The racing diffusion model
#'
#' A race between `n` independent diffusions, for choices with more than
#' two alternatives. Each accumulator is a Wiener process with its own
#' positive drift rate and unit diffusion, starting from a point drawn
#' uniformly on `(0, A)` and running to a common threshold `b`. The
#' first accumulator to reach the threshold is the response, and the
#' observed time is its arrival time plus a non-decision time.
#'
#' The model is [lba()]'s geometry with the ballistic assumption
#' removed: same start-point range, same threshold, but each accumulator
#' is noisy WITHIN a trial rather than drawing one rate at the start of
#' it. One accumulator's arrival time is therefore an inverse Gaussian
#' rather than the reciprocal of a normal, and the whole race stays
#' closed form.
#'
#' Use it rather than [wiener()] when a design has more than two
#' response alternatives, and rather than [lba()] when the story you
#' want for trial-to-trial spread is within-trial noise rather than a
#' rate drawn afresh each trial.
#'
#' @section The model:
#' Write the decision time as \eqn{t = rt - ndt}. One accumulator with
#' drift \eqn{v} travels a distance \eqn{b - z} from its start point
#' \eqn{z}, and \eqn{z} is uniform on \eqn{(0, A)}, so the DISTANCE is
#' uniform on \eqn{(k, k + A)} with \eqn{k = b - A}. For a fixed
#' distance \eqn{c} the first-passage time is inverse Gaussian with
#' density \eqn{c\,\phi((c - vt)/\sqrt{t})\,t^{-3/2}}, and averaging that
#' over the distance integrates in closed form. Put
#'
#' \deqn{x_1 = \frac{k - v t}{\sqrt t}, \qquad
#'       x_2 = \frac{k + A - v t}{\sqrt t},}
#'
#' the standardized positions of the two ends of the distance range at
#' time \eqn{t}. Then that accumulator's defective density and survival
#' are
#'
#' \deqn{f(t) = \frac{v\,(\Phi(x_2) - \Phi(x_1)) -
#'                    (\phi(x_2) - \phi(x_1))/\sqrt t}{A},}
#' \deqn{S(t) = \frac{\sqrt t\,(x_2\,\Delta\Phi + \Delta\phi) +
#'                    A\,\Phi(x_1) -
#'                    (\Delta E + \Delta\Phi)/(2 v)}{A},}
#'
#' where \eqn{\Delta} is the change of a quantity across the distance
#' range and \eqn{E(c) = e^{2 v c}\,\Phi(-(v t + c)/\sqrt t)} is the
#' reflected term of the inverse-Gaussian distribution function. A trial
#' on which accumulator \eqn{j} responded at time \eqn{rt} contributes
#'
#' \deqn{\log f_j(t) + \sum_{i \neq j} \log S_i(t),}
#'
#' the same race [lba()] runs, over a different single-accumulator law.
#'
#' @section Parameters:
#' \describe{
#'   \item{`v1`, ..., `vn`}{Drift rates, one per accumulator. LOG link,
#'     which is where this family parts company with [lba()]: an LBA
#'     drift is the MEAN of a normal drawn once per trial and may be
#'     negative, but a racing-diffusion drift is the rate itself, and an
#'     accumulator with a rate of zero or less never reaches the
#'     threshold at all. These are the primary parameters, so the main
#'     formula goes to all of them and each gets its own coefficients.
#'     Give one its own formula to move it alone.}
#'   \item{`A`}{Upper end of the start-point range, so start points are
#'     uniform on `(0, A)`. Log link.}
#'   \item{`k`}{Distance from the top of the start-point range to the
#'     threshold, so that `b = A + k`. Log link.}
#'   \item{`ndt`}{Non-decision time. Bounded link; see below.}
#' }
#'
#' The threshold is `A + k` rather than a free `b`, for the reason
#' [lba()] gives: a log link on `k` makes `b > A` structural, so a
#' threshold inside the start-point range, where a fraction of trials
#' would begin already finished, is not a state the optimizer can reach.
#'
#' `A` divides the density, so it cannot be zero. The log link keeps it
#' positive at every value of its linear predictor, and a race with no
#' start-point variability is a limit this family approaches rather than
#' a model it fits.
#'
#' @section Mapping to other software and to the paper:
#' The same model is written three ways. This family's spelling is
#' [lba()]'s, because the geometry is the same and moving between the
#' two families should cost one word.
#'
#' | this family | EMC2 | Tillman et al. (2020) |
#' |---|---|---|
#' | `v1`..`vn` | `v` | drift rate |
#' | `A` | `A` | start-point range |
#' | `k` | `B` | threshold less the range |
#' | `A + k` | `B + A` | threshold |
#' | `ndt` | `t0` | non-decision time |
#' | fixed at 1 | `s` | diffusion coefficient |
#'
#' EMC2's `B` is this family's `k` and NOT the threshold: its own
#' sampler draws the distance to travel as `B + runif(1, 0, A)`, so `B`
#' is the gap above the start-point range, exactly as `k` is here.
#' Reading it as the threshold moves every threshold by `A`.
#'
#' The diffusion coefficient is fixed at one, which is what identifies
#' the scale: multiplying `A`, `k` and every drift by one constant, and
#' the diffusion coefficient with them, leaves the distribution of
#' `(choice, rt)` unchanged. EMC2 divides `A`, `B` and `v` by its `s`
#' and so fixes the same quantity in the same place.
#'
#' @section Non-decision time:
#' The density is zero at and below `ndt`, so the likelihood has a hard
#' edge at `ndt = min(rt)`. As [wiener()] and [lba()] do, `ndt` gets a
#' logit scaled onto `(0, max_ndt)` rather than a log link, which makes
#' the constraint structural. `max_ndt` defaults to the smallest
#' observed response time, taken when the model frame is assembled. Pass
#' it explicitly to pin the bound, which matters if you will `predict()`
#' on new data whose minimum differs.
#'
#' @section The response, and why not `dec()`:
#' A trial is a `(choice, time)` pair. The time is the response and the
#' choice is per-row data, which reaches the family through `vint()`:
#'
#' ```
#' frm(bf(rt | vint(choice) ~ cond), family = rdm(3), data = dat)
#' ```
#'
#' `vint1` is the ACCUMULATOR INDEX: a whole number in `1..n` naming
#' which accumulator reached the threshold, counting from one. This is
#' [lba()]'s rule, including that a factor is not accepted, so recode it
#' with `as.integer(factor(choice))` and check that the level order
#' matches the accumulator numbering.
#'
#' `dec()` does not apply, for [lba()]'s reason: that term carries a 0/1
#' indicator naming one of two boundaries, and a race of `n`
#' accumulators needs `1..n`.
#'
#' @section Censoring and truncation:
#' A race has a survivor function in closed form, so `cens()` and
#' `trunc()` both work.
#'
#' The reason is that the race is over as soon as ANY accumulator
#' finishes, so the probability that the response time is past `t` is
#' the probability that none of them has:
#'
#' \deqn{S(t) = \prod_i S_i(t),}
#'
#' one factor per accumulator, each the start-point integral the
#' likelihood already forms for the losers of an observed trial. The
#' family declares that product as its log survivor function and its
#' complement as its distribution function, so nothing new is derived
#' for censoring: it is the same `lsurv` the density's loser terms use,
#' multiplied over all `n` accumulators instead of over `n - 1` of them.
#'
#' ```
#' frm(bf(rt | vint(choice) + cens(censored) ~ cond), family = rdm(3),
#'     data = dat)
#' ```
#'
#' A censored trial has no winner to report, because the race had not
#' finished when the clock ran out. `vint()` is still required, since a
#' declaration cannot be conditional on a censoring code, so give such a
#' row any accumulator index: the likelihood does not read it. The
#' distribution function is written as `-expm1(log S)`, which keeps its
#' digits where `1 - S` would lose them, and the log survivor goes to
#' `frmtmb` on the LOG scale, so a right-censored row stays exact past
#' the point where `log(1 - F)` is a constant with a zero gradient.
#'
#' @section Accuracy:
#' Every piece is written in the form that keeps its digits rather than
#' the form the paper prints, and the reasons are [lba()]'s reasons.
#'
#' The normal difference \eqn{\Phi(x_2) - \Phi(x_1)} goes through the
#' same log-space helper the linear ballistic accumulator uses, because
#' the subtractive spelling returns exactly zero once the lower tail
#' rounds to one. The density difference \eqn{\phi(x_2) - \phi(x_1)} and
#' the boundary-term difference \eqn{\Delta E} both go through
#' `ddm_expdiff()`, which anchors at the larger of the two logs:
#' written directly, the first underflows one term to zero while the
#' other is still finite, and the second overflows.
#'
#' The survival is written out rather than as `1 - F`. Measured against
#' a Gauss-Legendre reference that never subtracts, on a grid of 648
#' points, the form here holds 1.3e-12 relative down to a survival of
#' 1e-3, 6.4e-11 down to 1e-15 and 4.1e-8 down to 1e-250, and returns no
#' exact zeros. `1 - EMC2:::pWald()` returns exactly zero on 35 of the
#' 315 points of a comparable grid, the first at a survival near 1e-13,
#' which would send a loser's log contribution to `-Inf` on an ordinary
#' row.
#'
#' The density agrees with `EMC2:::dWald()` to 2.1e-4 relative and with
#' the same non-subtracting reference to 1.6e-13. The gap is the
#' subtractive normal difference in EMC2's form, not a disagreement
#' about the model.
#'
#' Below all of this the density underflows in double precision; the log
#' density is floored rather than returning `-Inf`, so the optimizer
#' sees a finite wall instead of a hole. Decision times at or below
#' zero, which `predict()` on faster new data can reach, are floored the
#' same way.
#'
#' @param n Number of accumulators, so the number of response
#'   alternatives. At least 2.
#' @param max_ndt Upper bound for the non-decision time, in the units of
#'   the response. `NULL`, the default, takes it from the data.
#'
#' @return A `frmtmb_family`.
#'
#' @references
#' Tillman, G., Van Zandt, T. and Logan, G. D. (2020). Sequential
#' sampling models without random between-trial variability: The racing
#' diffusion model of speeded decision making. *Psychonomic Bulletin and
#' Review*, 27(5), 911-936.
#'
#' @seealso [rdm_simulate()] to generate from the model, [lba()] for the
#'   ballistic race with the same geometry, and [wiener()] for the
#'   two-choice diffusion.
#'
#' @examples
#' set.seed(1)
#' dat <- rdm_simulate(400, v = c(3.0, 2.0, 1.2), A = 0.5, k = 0.5,
#'                     ndt = 0.2)
#' fit <- frm(bf(rt | vint(choice) ~ 1), family = rdm(3), data = dat)
#' fixef(fit)
#'
#' @export
rdm <- function(n, max_ndt = NULL) {
  if (missing(n) || !is.numeric(n) || length(n) != 1L || is.na(n) ||
      n != round(n) || n < 2) {
    stop("rdm(): `n` is the number of accumulators, one whole number ",
         "of 2 or more, e.g. rdm(3) for a three-alternative choice.",
         call. = FALSE)
  }
  n <- as.integer(n)
  if (!is.null(max_ndt)) {
    if (!is.numeric(max_ndt) || length(max_ndt) != 1L ||
        !is.finite(max_ndt) || max_ndt <= 0) {
      stop("rdm(): `max_ndt` bounds the non-decision time and must be ",
           "one positive finite number, or NULL to read it off the ",
           "response.", call. = FALSE)
    }
  }

  vp <- paste0("v", seq_len(n))
  dpn <- c(vp, "A", "k", "ndt")
  lk <- rep(list("log"), length(dpn))
  names(lk) <- dpn

  fam <- frmtmb::custom_family(
    "rdm",
    dpars = dpn,
    links = lk,
    lpdf = function(y, dpars, aterms) {
      lba_race_lpdf(y - dpars[["ndt"]], aterms[["vint1"]],
                    rdm_law, rdm_pars(dpars, vp))
    },
    valid_y = function(y, aterms) rdm_check_response(y, aterms, n),
    lccdf = function(q, dpars, aterms) rdm_lccdf(q, dpars, vp),
    # -expm1 rather than 1 - exp: the survivor of a race is within a
    # rounding of one for any `q` a fit spends time near, and the
    # subtractive form returns exactly zero there.
    lcdf = function(q, dpars, aterms) {
      ddm_floor(-expm1(rdm_lccdf(q, dpars, vp)), 1e-300)
    },
    family_finalize = function(fam, y, aterms) {
      ddm_ndt_finalize(fam, y, max_ndt, "rdm")
    },
    required_aterms = "vint1",
    init_dpars = rdm_inits(vp),
    type = "continuous",
    # A refusing mean rather than no mean at all. With the slot empty,
    # frmtmb falls back to "the first primary dpar on the response
    # scale is the mean", and for this family that is a DRIFT RATE:
    # predict(type = "response") returned 3.51 for a data set whose
    # response times average 0.36, and residuals(type = "response")
    # came back empty rather than refusing. A number that is not the
    # quantity it is labelled with is worse than an error.
    post = list(mean_fn = function(dpars, aterms) {
      stop("rdm: this family has no mean response time to report. The ",
           "mean of the winning accumulator's arrival is an ",
           "expectation over the minimum of several inverse-Gaussian ",
           "first passages and has no closed form, so fitted(), ",
           "predict(type = \"response\") and residuals(type = ",
           "\"response\") are all unavailable. predict(type = \"link\") ",
           "gives the drift rates, thresholds and non-decision time ",
           "the fit actually estimated.", call. = FALSE)
    }),
    sim = function(dpars, aterms, n_) {
      rdm_sim_rt(dpars, aterms, n_, vp)
    },
    primary_dpars = vp)

  # Carried so that a reader of the fitted object can see how many
  # accumulators raced, and so that the post-fit helpers do not have to
  # be told again.
  fam[["rdm_n"]] <- n
  fam
}

#' The one accumulator law the race is written against.
#'
#' The same two-function interface `lba_law` has, so `lba_race_lpdf()`
#' runs this race without knowing which law it is racing. `p` carries one
#' accumulator's drift `v`, start-point range `A` and threshold gap `k`.
#'
#' @noRd
rdm_law <- list(
  ldens = function(t, p) {
    st <- sqrt(t)
    x1 <- (p$k - p$v * t) / st
    x2 <- (p$k + p$A - p$v * t) / st
    dPhi <- lba_phidiff(x2, x1)
    dphi <- ddm_expdiff(ddm_lphi(x2), ddm_lphi(x1))
    f <- (p$v * dPhi - dphi / st) / p$A
    # Floored rather than clamped at zero, for the reason lba_law gives:
    # the race multiplies this by a 0/1 mask, and 0 * -Inf is NaN, which
    # would take the whole tape with it on the first underflowed row.
    log(ddm_floor(f, 1e-300))
  },
  lsurv = function(t, p) {
    st <- sqrt(t)
    x1 <- (p$k - p$v * t) / st
    x2 <- (p$k + p$A - p$v * t) / st
    dPhi <- lba_phidiff(x2, x1)
    dphi <- ddm_expdiff(ddm_lphi(x2), ddm_lphi(x1))
    # The integral of Phi over the distance range, regrouped so that the
    # only subtraction left is the one lba_phidiff() already protects.
    # The textbook form is a difference of two values of
    # x Phi(x) + phi(x), and once the whole range has moved past the
    # threshold those two collapse to the same tiny number.
    i1 <- st * (x2 * dPhi + dphi) + p$A * RTMB::pnorm(x1)
    # The reflected term of the inverse-Gaussian distribution function,
    # at each end of the distance range, in logs. It is part of a
    # probability and so is bounded above by one, which is why its log
    # is at most zero and neither exponential below can overflow.
    le1 <- 2 * p$v * p$k +
      RTMB::pnorm(-(p$v * t + p$k) / st, log.p = TRUE)
    le2 <- 2 * p$v * (p$k + p$A) +
      RTMB::pnorm(-(p$v * t + p$k + p$A) / st, log.p = TRUE)
    # The 1 / (2 v) is a removable singularity, not a division by
    # something that reaches zero: the bracket is the change across the
    # range of a quantity that is identically one at zero drift, so it
    # vanishes linearly in v and the ratio stays finite. Dividing rather
    # than expanding is what keeps the result EXACT away from zero,
    # which is where every fit lives, and the log link on the drifts is
    # what keeps v off zero itself.
    s <- (i1 - (ddm_expdiff(le2, le1) + dPhi) / (2 * p$v)) / p$A
    log(ddm_floor(s, 1e-300))
  })

#' Split the distributional parameters into one list per accumulator.
#'
#' @noRd
rdm_pars <- function(dpars, vp) {
  out <- vector("list", length(vp))
  for (j in seq_along(vp)) {
    out[[j]] <- list(v = dpars[[vp[j]]], A = dpars[["A"]],
                     k = dpars[["k"]])
  }
  out
}

#' Log probability that the race is still running at time `q`.
#'
#' The whole censoring seam, and it is three lines because the race
#' already computes what it needs: the trial is unfinished exactly when
#' every accumulator is, the accumulators are independent, and one
#' accumulator's survival is the start-point integral `rdm_law$lsurv()`
#' forms for every loser of every observed trial. So the log survivor of
#' the race is the sum of the same logs over ALL `n` accumulators rather
#' than over the `n - 1` that lost.
#'
#' `q` is a response time and the decision time is floored the way
#' `lba_race_lpdf()` floors it, for the same reason: at or below the
#' non-decision time the standardized positions change sign and the law
#' is meaningless rather than small. The floored value returns a log
#' survivor of zero, which is the right answer there - nothing can have
#' finished before the non-decision time.
#'
#' @noRd
rdm_lccdf <- function(q, dpars, vp) {
  t <- ddm_floor(q - dpars[["ndt"]], 1e-12)
  accs <- rdm_pars(dpars, vp)
  ls <- 0 * t
  for (j in seq_along(accs)) ls <- ls + rdm_law$lsurv(t, accs[[j]])
  ls
}

#' Response and choice validation.
#'
#' @noRd
rdm_check_response <- function(y, aterms, n) {
  ddm_refuse_dec("rdm", n, aterms)
  if (any(!is.finite(y)) || any(y <= 0)) {
    stop("rdm: the response must be a strictly positive, finite ",
         "response time. A time of zero or less leaves no decision ",
         "time for any non-decision time at all.", call. = FALSE)
  }
  ch <- aterms[["vint1"]]
  if (any(!is.finite(ch)) || any(ch != round(ch)) ||
      any(ch < 1) || any(ch > n)) {
    bad <- unique(ch[!is.finite(ch) | ch != round(ch) | ch < 1 | ch > n])
    stop("rdm(", n, "): the vint() choice indicator names which ",
         "accumulator reached the threshold and must be a whole number ",
         "from 1 to ", n, ". Saw ",
         paste(utils::head(sort(bad), 5), collapse = ", "),
         ". A factor is not accepted; recode it with ",
         "as.integer(factor(choice)) and check that the level order ",
         "matches the accumulator numbering.", call. = FALSE)
  }
  invisible(NULL)
}

#' Starting values.
#'
#' The drifts start apart rather than equal, for the reason
#' `lba_inits()` gives: all-equal drifts are a saddle for a race, since
#' no gradient then distinguishes the accumulators.
#'
#' These are RESPONSE-scale values, which is what `init_dpars` takes:
#' frmtmb applies the link itself. Writing them on the log scale instead
#' makes every one of them negative, so the log link returns `NaN` and
#' the start is silently discarded.
#'
#' The spacing is multiplicative rather than subtractive, so that it
#' stays positive however many accumulators race; a subtractive step
#' walks a late accumulator's start through zero, which the log link
#' cannot represent.
#'
#' @noRd
rdm_inits <- function(vp) {
  out <- list(A = function(y, aterms) 0.3,
              k = function(y, aterms) 0.5,
              ndt = function(y, aterms) 0.5 * min(y))
  for (j in seq_along(vp)) {
    out[[vp[j]]] <- local({
      jj <- j
      function(y, aterms) 2 * 0.9^(jj - 1)
    })
  }
  out[c(vp, "A", "k", "ndt")]
}

#' Draw response times conditional on each row's observed choice.
#'
#' `simulate()` replaces the response only and the choice is data, so the
#' joint draw is conditioned on it by rejection. See [rdm_simulate()] for
#' the joint draw, which is what a user usually wants.
#'
#' @noRd
rdm_sim_rt <- function(dpars, aterms, n_, vp) {
  ch <- rep(as.numeric(aterms[["vint1"]]), length.out = n_)
  A <- rep(as.numeric(dpars[["A"]]), length.out = n_)
  k <- rep(as.numeric(dpars[["k"]]), length.out = n_)
  ndt <- rep(as.numeric(dpars[["ndt"]]), length.out = n_)
  V <- vapply(vp, function(p) rep(as.numeric(dpars[[p]]), length.out = n_),
              numeric(n_))
  dim(V) <- c(n_, length(vp))
  out <- rep(NA_real_, n_)
  todo <- seq_len(n_)
  for (round in 1:200) {
    if (!length(todo)) break
    dr <- rdm_race_draw(A[todo], k[todo], V[todo, , drop = FALSE])
    hit <- dr$choice == ch[todo]
    out[todo[hit]] <- dr$time[hit] + ndt[todo[hit]]
    todo <- todo[!hit]
  }
  if (length(todo)) {
    stop("rdm: could not draw a response time for ", length(todo),
         " row(s) whose observed choice the fitted parameters almost ",
         "never produce. simulate() conditions each draw on that row's ",
         "vint() choice, so a choice the model gives a vanishing ",
         "probability has no draw to give. Use rdm_simulate() for an ",
         "unconditional draw of choice and time together.",
         call. = FALSE)
  }
  out
}

#' One unconditional pass of the race, from the generative process.
#'
#' Draws a start point, hence a distance to travel, then that distance's
#' first-passage time from the inverse Gaussian. Nothing here evaluates
#' the density: the simulator and the likelihood are two statements of
#' the same model, and a simulator written from the density could not
#' catch an error in it.
#'
#' @noRd
rdm_race_draw <- function(A, k, V) {
  m <- nrow(V)
  nacc <- ncol(V)
  tim <- matrix(0, m, nacc)
  for (j in seq_len(nacc)) {
    tim[, j] <- rdm_ig(k + stats::runif(m, 0, A), V[, j])
  }
  w <- max.col(-tim, ties.method = "first")
  list(choice = w, time = tim[cbind(seq_len(m), w)])
}

#' First-passage times of a unit-diffusion Wiener process to a distance.
#'
#' The inverse Gaussian with mean `dist / v` and shape `dist^2`, drawn by
#' Michael, Schucany and Haas's transformation: one chi-square with a
#' single degree of freedom gives a root of the defining quadratic, and
#' one uniform chooses between that root and its reflection.
#'
#' A drift at or below zero never arrives, and `Inf` is the honest time
#' for it: it loses every race it is in.
#'
#' @noRd
rdm_ig <- function(dist, v) {
  n <- length(dist)
  mu <- dist / v
  lam <- dist * dist
  y <- stats::rnorm(n)^2
  x <- mu + mu * mu * y / (2 * lam) -
    (mu / (2 * lam)) * sqrt(4 * mu * lam * y + mu * mu * y * y)
  u <- stats::runif(n)
  out <- ifelse(u <= mu / (mu + x), x, mu * mu / x)
  ifelse(v > 0, out, Inf)
}

#' Simulate from a racing diffusion model
#'
#' Draws from the generative process directly: a uniform start point per
#' accumulator, then that accumulator's inverse-Gaussian first-passage
#' time, and whichever arrives first. This is the joint draw of choice
#' and time, which is what the model produces and what a fit needs;
#' `simulate()` on a fitted object instead redraws only the time, holding
#' each row's observed choice.
#'
#' @param n Number of trials.
#' @param v Drift rates, one per accumulator, all positive. Its length
#'   sets the number of accumulators. May be a matrix with one row per
#'   trial, for a drift that varies with a covariate.
#' @param A Upper end of the start-point range.
#' @param k Distance from the top of the start-point range to the
#'   threshold, so the threshold is `A + k`.
#' @param ndt Non-decision time.
#'
#' @return A data frame with `choice` (the accumulator that reached the
#'   threshold, from 1) and `rt`.
#'
#' @examples
#' set.seed(1)
#' dat <- rdm_simulate(500, v = c(3.0, 2.0, 1.2), A = 0.5, k = 0.5,
#'                     ndt = 0.2)
#' table(dat$choice)
#' tapply(dat$rt, dat$choice, mean)
#'
#' @export
rdm_simulate <- function(n, v, A = 0.5, k = 0.5, ndt = 0.2) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("rdm_simulate(): `n` is the number of trials to draw, one ",
         "whole number of 1 or more.", call. = FALSE)
  }
  n <- as.integer(n)
  V <- if (is.matrix(v)) v else matrix(as.numeric(v), n, length(v),
                                       byrow = TRUE)
  if (nrow(V) != n) {
    stop("rdm_simulate(): a matrix `v` gives one row of drift rates per ",
         "trial, so it needs ", n, " rows, not ", nrow(V), ".",
         call. = FALSE)
  }
  nacc <- ncol(V)
  if (nacc < 2) {
    stop("rdm_simulate(): `v` needs a drift rate for each of at least ",
         "two accumulators; a race of one has nothing to lose to.",
         call. = FALSE)
  }
  if (any(!is.finite(V)) || any(V <= 0)) {
    stop("rdm_simulate(): every drift rate must be positive and finite. ",
         "A racing-diffusion drift is the rate itself, not the mean of ",
         "one, and an accumulator with a rate of zero or less never ",
         "reaches the threshold.", call. = FALSE)
  }
  A <- rep(as.numeric(A), length.out = n)
  k <- rep(as.numeric(k), length.out = n)
  if (any(A <= 0) || any(k <= 0)) {
    stop("rdm_simulate(): `A` and `k` must both be positive, so that ",
         "the threshold A + k lies above the start-point range and no ",
         "trial begins already finished.", call. = FALSE)
  }
  dr <- rdm_race_draw(A, k, V)
  data.frame(choice = dr$choice,
             rt = dr$time + rep(as.numeric(ndt), length.out = n))
}
