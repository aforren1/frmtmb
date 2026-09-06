# The probability that a Wiener process has NOT crossed the upper
# boundary by a given time. That one quantity is the whole no-go branch
# of wiener_gng(), and it is the reason this file exists: wiener()
# declares no lcdf, and its compatibility table says so in as many
# words, because the first-passage distribution function is a third
# series with its own truncation problem.
#
# Two series again, and blended the same way the density blends its two,
# for the same reason: the normalized time u = t / a^2 decides which one
# converges, u depends on a PARAMETER, and an RTMB tape records no
# comparisons. Both are evaluated at a fixed generous truncation and
# combined with a logistic weight in log(u).
#
# What is NOT the same as the density is the centre of the weight, and
# the difference is deliberate. The large-time route computes the no-go
# probability DIRECTLY, as the gambler's-ruin probability plus a
# correction; the small-time route computes it as 1 - F_upper. The
# second form cancels wherever a no-go trial is surprising, which is
# exactly where the likelihood is being pushed hardest, so the weight is
# centred far lower than the density's and hands over to the large-time
# route as soon as that route is accurate.
#
# HOW FAR LOWER is the whole question, and the first answer here was
# wrong. The centre shipped at 0.06, chosen on a grid whose relative
# start point stopped at 0.75. Widened to include w = 0.9, the blend
# reaches 7.75e-05 relative at t = 2.5, v = 5, a = 4, w = 0.90, where
# the true no-go probability is 8.2e-16. It is not the series and it is
# not the reference: at that row the large-time route ALONE is 6.7e-15,
# and the small-time route is 12.1 relative because it forms
# 1 - F_upper and there is nothing left to subtract from. At u = 0.156
# the 0.06 weight still gives that route a share of about 3e-05, and a
# tiny share of a hopeless number is what the 7.75e-05 is.
#
# The lesson is that a smooth blend has no safe side unless BOTH
# branches stay merely inaccurate outside their regime. This one has a
# branch that goes catastrophically wrong rather than gently wrong, so
# the weight has to reach saturation before that branch does.

# ------------------------------------------------------- tuning constants
#
# ddm_cdf_ks: the small-time image sum runs j = -ddm_cdf_ks .. ddm_cdf_ks.
#
# It was 12 at 0.3.0 and is 4, and the change is free rather than a
# trade. A term at index j carries exp(-2 j (w + j) / u), and the blend
# gives the small route a non-zero weight only below u = 0.197, where
# tanh has not yet saturated; at that u the j = 2 term is already
# exp(-2 * 2 * 2.5 / 0.197), which is 9e-23, and every larger j is
# smaller still. Measured over 840 rows spanning t in {0.05 .. 15}, v in
# {-2 .. 5}, a in {0.8 .. 4} and w in {0.25 .. 0.9}, every truncation
# from 2 to 12 gives BIT-IDENTICAL values AND bit-identical gradients:
# 840 of 840 rows equal, worst absolute log difference 0, worst gradient
# difference 0. Four is chosen rather than two so that the margin over
# the terms that do matter, |j| = 1, is not one term wide.
#
# The reason it is worth doing is the go/no-go variability integral,
# which evaluates this function on a three-dimensional node grid. The
# image sum is where the pnorm calls are, so 25 terms against 9 is most
# of the cost of every node.
#
# ddm_cdf_kl: the large-time series runs k = 1 .. ddm_cdf_kl. NOT
# reducible the same way: its terms decay like exp(-k^2 pi^2 u / 2), and
# the blend hands over to it from u = 0.01, where k = 30 is only
# exp(-44).
# ddm_cdf_u0, ddm_cdf_us: centre and scale of the logistic weight in
# log(u).
#
# u0 is MEASURED, on 1200 rows spanning t in {0.05 .. 15}, v in
# {-2, -0.5, 0.5, 1, 2, 5}, a in {0.8, 1.4, 2.5, 4} and w in
# {0.25, 0.45, 0.5, 0.75, 0.9}, against a 260-bit reference:
#
#   u0     us     max rel    rows >1e-12   >1e-9
#   0.06   0.12   7.75e-05        14         5     <- shipped at 0.3.0
#   0.04   0.12   9.00e-08         7         1
#   0.03   0.12   7.45e-10         3         0
#   0.02   0.12   2.17e-12         1         0     <- here
#   0.015  0.12   2.17e-12         2         0
#   0.01   0.12   4.57e-11         4         0
#   0.005  0.12   3.24e-06        29         4
#
# Lowering it further is not free, which is why 0.02 and not 0.005: the
# two routes cross near u = 0.025, and below the crossing the SMALL
# route is the accurate one, so an over-low centre hands over too early
# and the error climbs again. 0.02 sits just under the crossing, which
# is the right side to be on because the small route's degradation
# above it is far steeper than the large route's below it.
ddm_cdf_ks <- 4L
ddm_cdf_kl <- 30L
ddm_cdf_u0 <- 0.02
ddm_cdf_us <- 0.12

# The largest exponent any term is allowed to reach before it is
# exponentiated. Not an accuracy knob: no term of a probability can
# legitimately be of order e^500, so this only ever fires where the
# series is outside its own regime and its weight is zero. It is here
# because an Inf on the tape poisons the reverse pass even where its
# weight is zero, which is the same reason the density uses tanh rather
# than the exp spelling of the same curve.
ddm_cdf_cap <- 500

#' Log of `sinh(x) / x`.
#'
#' The gambler's-ruin probability is a ratio of two exponential
#' differences that is 0/0 at zero drift. Factoring each difference as
#' `2 exp(mid) sinh(half)` moves the whole singularity into
#' `sinh(x) / x`, which is smooth and equal to one there.
#'
#' Written as a LOG, and through `log1p`, so that it neither overflows
#' nor loses the small-argument limit. `sinh(x) / x` itself overflows
#' once `|x|` passes about 710, which a boundary separation times a
#' drift can reach on a wild optimizer step; the log cannot.
#'
#' The regularizer is 1e-20 rather than `ddm_pos()`'s 1e-40, and the
#' difference matters. `ddm_pos()` would give `r` around 1e-20 at zero,
#' where `exp(-2 r)` rounds to exactly one and `log1p(-1)` is `-Inf`.
#' At 1e-20 the smallest `r` is 1e-10, `exp(-2 r)` is safely below one,
#' and the perturbation is invisible: `log sinh(x)/x` is quadratic in
#' `x` near zero, so an absolute error of 1e-10 in `r` moves the answer
#' by about 1e-21.
#'
#' @noRd
ddm_lsinhc <- function(x) {
  r <- sqrt(x * x + 1e-20)
  r + log1p(-exp(-2 * r)) - log(2 * r)
}

#' Log of the probability of ever being absorbed at the LOWER boundary.
#'
#' The classical gambler's-ruin result,
#' `(exp(-2 v a) - exp(-2 v a w)) / (exp(-2 v a) - 1)`, in the one
#' spelling that survives both of its failure modes: the ratio is 0/0 at
#' zero drift, and both exponentials overflow for a large enough drift
#' times boundary separation. Factoring through `sinh` removes the
#' first, and taking the log removes the second.
#'
#' `log1p(-w)` rather than `log(1 - w)` because `w` reaches the
#' neighbourhood of one whenever the start point sits near the upper
#' boundary, which is an ordinary bias and not an extreme.
#'
#' @noRd
ddm_llower_prob <- function(v, a, w) {
  -v * a * w + log1p(-w) +
    ddm_lsinhc(v * a * (1 - w)) - ddm_lsinhc(v * a)
}

#' Defective distribution function of the LOWER boundary, small time.
#'
#' The image series of Blurton, Kesselmeier and Vollrath (2012), reached
#' by integrating the small-time density term by term. Each term of that
#' density is an inverse-Gaussian density at a SIGNED distance
#' `c_j = a (w + 2 j)`, so each integrates to an inverse-Gaussian
#' distribution function
#'
#'   A(t, c) = Phi((v t - c)/sqrt(t)) + exp(2 v c) Phi(-(v t + c)/sqrt(t))
#'
#' The integration constant is what differs between the terms: `A` tends
#' to zero as `t` does when `c` is positive, and to `1 + exp(2 v c)` when
#' `c` is negative. Since `w` is in `(0, 1)`, the sign of `c_j` is
#' decided by `j` alone. That is what makes this tape-safe: the branch is
#' on the LOOP INDEX, which is structure, not on a parameter.
#'
#' Every term is one `exp()` of one summed exponent, never a product of
#' an exponential and a probability, because the two factors are of
#' opposite and enormous magnitude for `|j|` large and only their sum is
#' of a size a double can hold.
#'
#' @noRd
ddm_lower_cdf_small <- function(t, v, a, w, K = ddm_cdf_ks) {
  st <- sqrt(t)
  out <- 0 * t
  for (j in 0:K) {
    cj <- a * (w + 2 * j)
    out <- out +
      exp(ddm_smin(-2 * v * a * (w + j) +
                     RTMB::pnorm((v * t - cj) / st, log.p = TRUE),
                   ddm_cdf_cap)) +
      exp(ddm_smin(2 * v * a * j +
                     RTMB::pnorm(-(v * t + cj) / st, log.p = TRUE),
                   ddm_cdf_cap))
  }
  for (j in (-K):(-1)) {
    cj <- a * (w + 2 * j)
    out <- out -
      exp(ddm_smin(-2 * v * a * (w + j) +
                     RTMB::pnorm(-(v * t - cj) / st, log.p = TRUE),
                   ddm_cdf_cap)) -
      exp(ddm_smin(2 * v * a * j +
                     RTMB::pnorm((v * t + cj) / st, log.p = TRUE),
                   ddm_cdf_cap))
  }
  out
}

#' No-go probability, small-time route.
#'
#' The upper boundary is the lower one reflected, so the go boundary's
#' distribution function is `ddm_lower_cdf_small()` at `-v` and `1 - w`.
#' Subtracting from one is accurate here and only here: at a small
#' normalized time the process has had no chance to reach the go
#' boundary, so the quantity being subtracted is small and the
#' difference keeps its digits.
#'
#' @noRd
ddm_nogo_small <- function(t, v, a, w, K = ddm_cdf_ks) {
  1 - ddm_lower_cdf_small(t, -v, a, 1 - w, K)
}

#' No-go probability, large-time route.
#'
#' Directly, rather than as `1 - F_upper`. The eventual lower-boundary
#' probability plus one correction series:
#'
#'   nogo(t) = P_lower + 2 pi exp(v a (1 - w))
#'             sum_k (-1)^(k+1) k sin(k pi w) exp(-lambda_k t)
#'                    / (v^2 a^2 + k^2 pi^2)
#'
#' with `lambda_k = (v^2 a^2 + k^2 pi^2) / (2 a^2)`. It is the sum of the
#' large-time series for the lower-boundary distribution function and
#' the one for the probability of no absorption at all, which share
#' every factor and collapse into a single alternating sum.
#'
#' Why it is worth having a second route at all: this one never forms
#' `1 - F`, so it keeps full relative precision on a no-go probability of
#' 1e-13, where the subtraction has none. That regime is not exotic. It
#' is a no-go trial under a model that expects a response, which is
#' precisely the row a fit is being pulled by.
#'
#' @noRd
ddm_nogo_large <- function(t, v, a, w, K = ddm_cdf_kl) {
  va2 <- v * v * a * a
  drift <- v * a * (1 - w)
  s <- 0 * t
  for (k in 1:K) {
    den <- va2 + k * k * pi * pi
    s <- s + ((-1)^(k + 1) * k * sin(k * pi * w) / den) *
      exp(ddm_smin(drift - den * t / (2 * a * a), ddm_cdf_cap))
  }
  exp(ddm_llower_prob(v, a, w)) + 2 * pi * s
}

#' Log probability of no UPPER-boundary crossing by time `t`.
#'
#' `t` is decision time, so the non-decision time is already subtracted.
#' `v`, `a` and `w` are the drift, boundary separation and relative start
#' point of [wiener()], unreflected: the upper boundary is the go
#' boundary and the reflection happens inside.
#'
#' Both probabilities are floored before their logs are taken. Outside
#' its own regime a series sum is not merely inaccurate, it can be
#' negative, and `log()` of a negative is `NaN` while `0 * NaN` is `NaN`
#' too, so a saturated weight alone would not keep the wrong branch out
#' of the answer. This is `ddm_log_gs()`'s argument for `ddm_pos()`,
#' applied to a probability instead of a series sum.
#'
#' @noRd
ddm_nogo_lprob <- function(t, v, a, w) {
  u <- ddm_floor(t / (a * a), ddm_u_floor)
  lam <- 0.5 * (1 + tanh((log(u) - log(ddm_cdf_u0)) / ddm_cdf_us))
  ls <- log(ddm_floor(ddm_nogo_small(t, v, a, w), ddm_share_floor))
  ll <- log(ddm_floor(ddm_nogo_large(t, v, a, w), ddm_share_floor))
  (1 - lam) * ls + lam * ll
}

# --------------------------------------------- across-trial variability
#
# The go/no-go likelihood needs the no-go probability averaged over the
# same three per-trial distributions [wiener()] averages its density
# over, and NONE of the three has the closed form the density's drift
# integral has.
#
# The drift is the one that is different, and the reason is worth
# stating because it decides the interface. In the DENSITY the drift
# appears only as `exp(-v a w - v^2 t / 2)`, an exponential-quadratic,
# so a normal average completes the square and costs one square root.
# In the DISTRIBUTION FUNCTION the drift also enters the eigenvalues,
# through `1 / (v^2 a^2 + k^2 pi^2)` in every term of the large-time
# series and through the gambler's-ruin probability the series corrects,
# and neither of those is an exponential-quadratic. So `sv` needs a
# quadrature here where it needs none there, and `wiener_gng(nodes =)`
# carries an `sv` entry that `wiener(nodes =)` has no use for.
#
# The three integrands are not equally hard, and the measured node
# counts in `gng_check_nodes()` say so.

#' No-go probability with the across-trial variability integrated out.
#'
#' `td` is the deadline and `t0` the mean non-decision time, both on the
#' response's scale: unlike `ddm_nogo_lprob()`, which takes a decision
#' time, this one subtracts the non-decision time itself because it is
#' integrating over it.
#'
#' `st` shifts the DEADLINE rather than cutting a range. The density's
#' non-decision-time integral has its range cut at the response time,
#' because past that the decision time is negative and the density is
#' identically zero; the no-go probability has no such cut, because a
#' non-decision time past the deadline leaves the accumulator no time at
#' all and the probability of no response is then exactly one. The whole
#' uniform range is therefore averaged, with the decision window floored
#' at `ddm_gng_t_floor`, where the small-time route returns one to every
#' bit a double holds.
#'
#' The average is over PROBABILITIES, so it is a log-sum-exp, anchored
#' at the maximum over the nodes for the reason `ddm_lpdf_var()` gives:
#' every exponent is then at most zero and one of them is exactly zero,
#' whatever the parameters do.
#'
#' @noRd
ddm_nogo_lprob_var <- function(td, v, a, w, t0, sv, sz, st, nd) {
  nv <- length(nd$sv$x)
  nz <- length(nd$sz$x)
  nt <- length(nd$st$x)
  lps <- vector("list", nv * nz * nt)
  wts <- numeric(length(lps))
  m <- 0L
  for (j in seq_len(nt)) {
    tau <- t0 + st * (nd$st$x[[j]] - 0.5)
    tt <- ddm_floor(td - tau, ddm_gng_t_floor)
    for (i in seq_len(nz)) {
      om <- w + ddm_wclamp_dev(sz * (nd$sz$x[[i]] - 0.5), w)
      for (h in seq_len(nv)) {
        m <- m + 1L
        lps[[m]] <- ddm_nogo_lprob(tt, v + sv * nd$sv$x[[h]], a, om)
        wts[m] <- nd$st$w[[j]] * nd$sz$w[[i]] * nd$sv$w[[h]]
      }
    }
  }
  if (length(lps) == 1L) return(lps[[1L]])
  lref <- lps[[1L]]
  for (m in seq_along(lps)[-1L]) {
    d <- lref - lps[[m]]
    lref <- 0.5 * (lref + lps[[m]] + abs(d))
  }
  acc <- 0
  for (m in seq_along(lps)) acc <- acc + wts[m] * exp(lps[[m]] - lref)
  lref + log(acc)
}

#' Hold a relative start point strictly inside the boundaries.
#'
#' A uniform start point of width `sz` centred on a `bias` near one runs
#' past the upper boundary, and [wiener()] says so: the parameterization
#' keeps the range inside the boundaries by construction only at an
#' unbiased start. What [wiener()] does NOT say is that the two branches
#' fail differently out there, and the difference is the whole reason
#' this function exists.
#'
#' Measured, at `t = 1.25`, drift 1, separation 1.4:
#'
#'   w      no-go log prob     density log
#'   0.95   -4.441             -5.988
#'   1.00   -41.04             -40.81
#'   1.05   NaN                -6.057
#'   -0.05  +0.1412            -5.089
#'
#' The density returns a finite number, which is wrong but is a barrier
#' an optimizer climbs out of. The probability returns `NaN` above one,
#' because the gambler's-ruin term takes `log1p(-w)`, and a positive log
#' probability below zero, which is a probability above one. `NaN` is
#' not a barrier; it is the end of the tape, and one node of one row is
#' enough to take the whole fit with it.
#'
#' The clamped value is the correct limit FOR THE BRANCH IT IS APPLIED
#' TO, and that qualifier is the whole of the honesty here. A start
#' point at the upper boundary is a trial that has already finished at
#' the go boundary, so its no-go probability is zero and `log` of it
#' runs to minus infinity; one at the lower boundary has finished at the
#' other, so its no-go probability is one. Both are what this returns.
#'
#' What it is NOT is a repair of the pair. The go branch is not clamped,
#' deliberately, so past a boundary the two branches average different
#' start-point distributions and the go mass plus the no-go probability
#' no longer sums to one: measured, 0.970639 at `sz` = 0.5 with
#' `bias` = 0.85 and 0.843370 at `sz` = 0.9 with `bias` = 0.90. Up to
#' 16 percent of the mass is gone.
#'
#' That is tolerable only because of what it is compared against. The
#' alternative here was a `NaN`, which takes the tape with it; and the
#' region is strictly downhill, so the clamp acts as a barrier rather
#' than an attractor. Profiled on 2000 trials at `bias` = 0.85 with a
#' true `sz` of 0.20, the best point inside the boundaries beats the
#' best point outside by 189 log units and the surface is monotone
#' across the crossing. No optimizer ends up in the defective region;
#' it walks away from it. `?wiener_gng` says all of this to users.
#'
#' The go branch is NOT clamped, and deliberately: it is [wiener()]'s
#' own averaged density, shared bit for bit with that family, and
#' clamping it here would make the same model score differently under
#' two family names. The asymmetry is confined to a region neither
#' branch models.
#'
#' It clamps the DEVIATION rather than the start point, and that is not
#' a stylistic choice. `ddm_floor()` and `ddm_smin()` are arithmetic, so
#' clamping the value perturbs it by an ulp even where the range fits
#' entirely inside the boundaries: measured, `ddm_wclamp(0.5)` came back
#' as 0.4999999999999999. That ulp costs two exact identities this
#' family is otherwise entitled to - the no-variability path is no
#' longer bit-identical to the plain one, and a right-censored row no
#' longer scores exactly as a no-go row - for no accuracy at all.
#' Clamping the deviation is EXACTLY inert at a deviation of zero,
#' because the whole expression collapses to `lo + (-lo)`, which is
#' exactly zero in IEEE arithmetic whatever `lo` is.
#'
#' @noRd
ddm_wclamp_dev <- function(d, w) {
  ddm_floor(ddm_smin(d, 1 - ddm_w_edge - w), ddm_w_edge - w)
}

# How far inside the boundaries a start point is held. Far enough below
# one that log1p(-w) is an ordinary -20.7 rather than -Inf, and far
# enough above zero that the same is true of log(w).
ddm_w_edge <- 1e-9

# The floor the go/no-go decision window is held at. A deadline that has
# not yet cleared the non-decision time leaves no window in which a
# response could happen, so the no-go probability there is exactly one;
# the small-time route returns one to every bit a double holds at this
# floor, and the floor keeps the log off a negative argument.
ddm_gng_t_floor <- 1e-12
