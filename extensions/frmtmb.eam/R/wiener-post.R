# Post-fit quantities and simulation. Everything here runs on doubles,
# off the AD tape, so it may branch freely.

#' Probability that the accumulator reaches the UPPER boundary.
#'
#' The classical gambler's-ruin answer,
#' `(1 - exp(-2 v a w)) / (1 - exp(-2 v a))`, written with `expm1()` so
#' that it stays accurate as the drift goes to zero, where both
#' differences vanish and the naive spelling loses every digit.
#'
#' @noRd
ddm_p_upper <- function(v, a, w) {
  m <- 2 * v * a
  # A ratio of two expm1() calls, each accurate to a rounding error of
  # its own result, so the ratio is accurate everywhere except at
  # exactly zero drift, where it is 0/0 and the limit is the relative
  # start point. No tolerance band is needed here; the conditional
  # means below are the ones that cancel.
  out <- expm1(-m * w) / expm1(-m)
  out[m == 0] <- rep_len(w, length(m))[m == 0]
  out
}

# Below this value of |2 * v * a| the conditional means use their
# zero-drift limits.
#
# Two errors meet here. The direct formulas subtract two quantities of
# order m^2 to get one of order m^3, twice over, so they lose accuracy
# like eps / m^2. The limits are wrong by order m^2, NOT order m: the
# first-order term of the expansion vanishes identically, which is why
# the band can be as wide as this without costing anything. Balancing
# 0.02 m^2 against 50 eps / m^2 puts the crossover near 1e-3 and the
# worst case either side of it below 1e-8 relative. Measured in
# tests/testthat/test-moments.R against quadrature of RWiener's density.
ddm_drift_tol <- 1e-3

#' Mean decision time conditional on the boundary reached.
#'
#' Solves the two standard boundary-value problems for the diffusion:
#' `T(z)` is the mean exit time from start `z`, and `u(z)` is the mean
#' exit time weighted by the event "exits at the upper boundary", so
#' `u / P(upper)` and `(T - u) / P(lower)` are the two conditional
#' means. At zero drift both reduce to the polynomials
#' `(a^2 - z^2) / 3` and `z (2a - z) / 3`.
#'
#' @noRd
ddm_cond_mean_dt <- function(v, a, w, up) {
  n <- length(v)
  z <- w * a
  m <- 2 * v * a
  out <- numeric(n)
  small <- abs(m) < ddm_drift_tol

  if (any(small)) {
    i <- which(small)
    out[i] <- ifelse(up[i] == 1,
                     (a[i]^2 - z[i]^2) / 3,
                     z[i] * (2 * a[i] - z[i]) / 3)
  }
  if (any(!small)) {
    i <- which(!small)
    vv <- v[i]; aa <- a[i]; ww <- w[i]; zz <- z[i]; mm <- m[i]
    emz <- exp(-mm * ww)          # exp(-2 v z)
    ema <- exp(-mm)               # exp(-2 v a)
    D <- -expm1(-mm)              # 1 - exp(-2 v a)
    piup <- expm1(-mm * ww) / expm1(-mm)
    tmean <- aa * (piup - ww) / vv
    C1 <- aa * (1 + ema) / (vv * D * D)
    uz <- -(zz / (vv * D)) * (1 + emz) + C1 * (-expm1(-mm * ww))
    out[i] <- ifelse(up[i] == 1, uz / piup, (tmean - uz) / (1 - piup))
  }
  out
}

#' The decision indicator for a mean, or the refusal.
#'
#' Two helpers rather than one with the noun as an argument, and two
#' rather than four inline copies. Every post-fit quantity here is
#' conditional on the boundary a row ended at, so each of them has
#' something to say when the boundary is absent, but a mean and a draw
#' do not have the SAME thing to say. One template per thing, each
#' written once, so a reported message resolves to one line of source.
#'
#' @noRd
ddm_indicator_mean <- function(aterms) {
  up <- ddm_indicator(aterms)
  if (is.null(up)) {
    stop("wiener: the mean response time is conditional on the boundary ",
         "a trial ended at, which is missing. See ?wiener.", call. = FALSE)
  }
  up
}

#' The decision indicator for a draw, or the refusal.
#'
#' @noRd
ddm_indicator_sim <- function(aterms) {
  up <- ddm_indicator(aterms)
  if (is.null(up)) {
    stop("wiener: simulating a response time needs the boundary each ",
         "row ended at, which is missing. See ?wiener.", call. = FALSE)
  }
  up
}

#' The family's `post$mean_fn`: expected response time.
#'
#' The response is the response time and the boundary is data, so the
#' mean frmtmb wants for a row is the mean response time GIVEN that
#' row's boundary: the conditional mean decision time plus the
#' non-decision time.
#'
#' @noRd
ddm_mean_rt <- function(dpars, aterms) {
  up <- ddm_indicator_mean(aterms)
  n <- max(lengths(list(dpars[["mu"]], dpars[["bs"]],
           dpars[["ndt"]], dpars[["bias"]], up)))
  v <- rep_len(dpars[["mu"]], n); a <- rep_len(dpars[["bs"]], n)
  t0 <- rep_len(dpars[["ndt"]], n); w <- rep_len(dpars[["bias"]], n)
  u <- rep_len(up, n)
  t0 + ddm_cond_mean_dt(v, a, w, u)
}

#' The family's `sim` slot: a response time for each row, given that
#' row's boundary.
#'
#' Exact inverse-transform sampling when RWiener is available.
#' `RWiener::qwiener()` inverts the DEFECTIVE distribution function, the
#' one that integrates to the boundary probability rather than to one,
#' so scaling a uniform draw by that probability samples the conditional
#' distribution directly. No rejection: a boundary the fitted parameters
#' almost never reach costs the same as a common one.
#'
#' @noRd
ddm_sim_rt <- function(dpars, aterms, n) {
  up <- ddm_indicator_sim(aterms)
  v <- rep_len(dpars[["mu"]], n); a <- rep_len(dpars[["bs"]], n)
  t0 <- rep_len(dpars[["ndt"]], n); w <- rep_len(dpars[["bias"]], n)
  u <- rep_len(up, n)

  if (!requireNamespace("RWiener", quietly = TRUE)) {
    return(ddm_sim_euler(v, a, w, t0, u, n))
  }
  pb <- ddm_p_upper(v, a, w)
  pb <- ifelse(u == 1, pb, 1 - pb)
  # qwiener returns Inf at the very top of the defective range, where
  # its bracketing search runs out; the clamp keeps the draw inside it
  p <- stats::runif(n) * pb * (1 - 1e-9)
  out <- numeric(n)
  for (i in seq_len(n)) {
    out[i] <- RWiener::qwiener(p[i], a[i], t0[i], w[i], v[i],
                               resp = if (u[i] == 1) "upper" else "lower")
  }
  bad <- !is.finite(out)
  if (any(bad)) {
    out[bad] <- ddm_sim_euler(v[bad], a[bad], w[bad], t0[bad], u[bad],
                              sum(bad))
  }
  out
}

#' Forward simulation of the diffusion, used when RWiener is absent.
#'
#' Euler-Maruyama on a fixed step with rejection to match the requested
#' boundary. Biased at order sqrt(dt) because a discretized path can
#' cross a boundary and come back between two steps, so it is the
#' fallback and not the default.
#'
#' @noRd
ddm_sim_euler <- function(v, a, w, t0, up, n, dt = 1e-4, tmax = 10,
                          max_pass = 50L) {
  out <- rep(NA_real_, n)
  todo <- seq_len(n)
  pass <- 0L
  while (length(todo)) {
    pass <- pass + 1L
    if (pass > max_pass) {
      stop("wiener: the fallback simulator could not produce a draw at ",
           "the requested boundary for ", length(todo), " of ", n,
           " rows in ", max_pass, " passes. Install RWiener for exact ",
           "conditional draws.", call. = FALSE)
    }
    k <- length(todo)
    x <- a[todo] * w[todo]
    live <- rep(TRUE, k)
    tt <- numeric(k); hitup <- integer(k)
    step <- 0L
    while (any(live) && step * dt < tmax) {
      step <- step + 1L
      x[live] <- x[live] + v[todo][live] * dt +
        stats::rnorm(sum(live), 0, sqrt(dt))
      hi <- live & x >= a[todo]
      lo <- live & x <= 0
      tt[hi | lo] <- step * dt
      hitup[hi] <- 1L
      live <- live & !hi & !lo
    }
    ok <- !live & hitup == up[todo]
    out[todo[ok]] <- tt[ok] + t0[todo][ok]
    todo <- todo[!ok]
  }
  out
}

#' Simulate a drift-diffusion data set
#'
#' Draws response times and boundary choices jointly from the Wiener
#' first-passage process, which is what an experiment produces: the
#' boundary is an outcome, not a design variable. The family itself
#' treats the boundary as data and simulates the response time
#' conditional on it, so this is the function to use when you want a
#' whole data set rather than a posterior predictive draw.
#'
#' Across-trial variability is drawn first and then handed to the same
#' per-trial draw, because the variability parameters do not change the
#' process: they say that each trial runs the ordinary diffusion with
#' its own drift rate, start point and non-decision time. That makes the
#' simulator an independent statement of the model from the density,
#' which is what a recovery study needs it to be.
#'
#' @param n Number of trials.
#' @param mu,bs,ndt,bias Drift rate, boundary separation, non-decision
#'   time and relative start point, on the response scale. Each may be a
#'   single value or a vector of length `n`, which is how a
#'   two-condition design is built.
#' @param sv,sz,st Across-trial variability: the standard deviation of a
#'   normal drift rate, and the widths of a uniform relative start point
#'   and a uniform non-decision time. Zero, the default for each, is the
#'   plain Wiener process.
#'
#' @return A data frame with `rt`, the response time, and `upper`, 1 for
#'   a response at the upper boundary and 0 for the lower one, in the
#'   coding `dec()` and `vint()` both expect.
#'
#' @examples
#' set.seed(1)
#' cond <- rep(c(0, 1), each = 100)
#' dat <- ddm_simulate(200, mu = 0.2 + 1.1 * cond, bs = 1.4,
#'                     ndt = 0.3, bias = 0.5)
#' dat$cond <- factor(cond)
#' str(dat)
#'
#' # Ratcliff's full model, at values the literature uses
#' full <- ddm_simulate(200, mu = 1.2, bs = 1.5, ndt = 0.3,
#'                      sv = 1.0, sz = 0.2, st = 0.1)
#'
#' @export
ddm_simulate <- function(n, mu, bs, ndt, bias = 0.5,
                         sv = 0, sz = 0, st = 0) {
  v <- rep_len(mu, n); a <- rep_len(bs, n)
  t0 <- rep_len(ndt, n); w <- rep_len(bias, n)
  if (any(a <= 0) || any(t0 < 0) || any(w <= 0) || any(w >= 1)) {
    stop("ddm_simulate(): need bs > 0, ndt >= 0 and bias strictly ",
         "inside (0, 1).", call. = FALSE)
  }
  sv <- rep_len(sv, n); sz <- rep_len(sz, n); st <- rep_len(st, n)
  if (any(sv < 0) || any(sz < 0) || any(st < 0)) {
    stop("ddm_simulate(): the across-trial variability widths sv, sz ",
         "and st are not negative.", call. = FALSE)
  }
  if (any(sv > 0)) v <- stats::rnorm(n, v, sv)
  if (any(sz > 0)) w <- stats::runif(n, w - sz / 2, w + sz / 2)
  if (any(st > 0)) t0 <- stats::runif(n, t0 - st / 2, t0 + st / 2)
  if (any(w <= 0) || any(w >= 1)) {
    stop("ddm_simulate(): sz pushed the relative start point outside ",
         "(0, 1); the uniform range must fit between the boundaries.",
         call. = FALSE)
  }
  if (any(t0 < 0)) {
    stop("ddm_simulate(): st pushed the non-decision time below zero; ",
         "the uniform range needs st / 2 no larger than ndt.",
         call. = FALSE)
  }
  if (requireNamespace("RWiener", quietly = TRUE)) {
    q <- numeric(n); r <- integer(n)
    for (i in seq_len(n)) {
      d <- RWiener::rwiener(1, a[i], t0[i], w[i], v[i])
      q[i] <- d$q
      r[i] <- as.integer(d$resp == "upper")
    }
    return(data.frame(rt = q, upper = r))
  }
  # joint draw from the fallback: simulate forward and keep the
  # boundary the path actually reached
  up <- stats::rbinom(n, 1, ddm_p_upper(v, a, w))
  data.frame(rt = ddm_sim_euler(v, a, w, t0, up, n), upper = up)
}

# ------------------------------------- post-fit under across-trial variability
#
# Everything below runs on doubles, off the tape, so it may branch and
# it may use rules the density cannot. The point of writing it at all is
# that the plain closed forms are the WRONG answers once the parameters
# vary between trials, and a post-fit method that quietly returns a
# number for the wrong model is the one failure this package works
# hardest to avoid.

#' Mean response time conditional on the boundary, with variability.
#'
#' Conditioning is what makes this more than an average. The trials that
#' reached a given boundary are not a fair sample of the drift rates and
#' start points: a drift toward that boundary reaches it more often. So
#' the conditional mean is a RATIO,
#'
#'   sum w P(boundary | nu, zeta) (t0 + T(nu, zeta))
#'   -------------------------------------------------
#'   sum w P(boundary | nu, zeta)
#'
#' with the boundary probability and the conditional mean decision time
#' both in closed form from `ddm_p_upper()` and `ddm_cond_mean_dt()`.
#'
#' The non-decision time drops out of the ratio: it is added to every
#' trial regardless of which boundary was reached, so its distribution
#' contributes its mean and nothing else. Only the drift and the start
#' point need nodes, and the drift's are Gauss-Hermite because its
#' distribution is normal.
#'
#' @noRd
ddm_mean_rt_var <- function(dpars, aterms, nd, gh) {
  up <- ddm_indicator_mean(aterms)
  n <- max(lengths(list(dpars[["mu"]], dpars[["bs"]], dpars[["ndt"]],
                        dpars[["bias"]], up)))
  g <- function(nm, default) {
    if (is.null(dpars[[nm]])) rep_len(default, n) else rep_len(dpars[[nm]], n)
  }
  v <- g("mu", 0); a <- g("bs", 1); t0 <- g("ndt", 0); w <- g("bias", 0.5)
  sv <- g("sv", 0); sz <- g("sz", 0)
  u <- rep_len(up, n)

  num <- numeric(n)
  den <- numeric(n)
  for (k in seq_along(gh$x)) {
    nu <- v + sv * gh$x[[k]]
    for (i in seq_along(nd$sz$x)) {
      om <- w + sz * (nd$sz$x[[i]] - 0.5)
      pu <- ddm_p_upper(nu, a, om)
      p <- ifelse(u == 1, pu, 1 - pu)
      wt <- gh$w[[k]] * nd$sz$w[[i]]
      num <- num + wt * p * ddm_cond_mean_dt(nu, a, om, u)
      den <- den + wt * p
    }
  }
  t0 + num / den
}

#' A conditional draw with the trial's own parameters.
#'
#' Drawing the per-trial parameters and then drawing a response time at
#' the requested boundary would be wrong for the same reason the mean is
#' a ratio: conditioning on the boundary reweights which drift rates and
#' start points the trial could have had. The draw is therefore accepted
#' with probability equal to the boundary probability it implies, which
#' samples that reweighting exactly, and only then is the response time
#' drawn from the plain conditional density.
#'
#' @noRd
ddm_sim_rt_var <- function(dpars, aterms, n, nd) {
  up <- ddm_indicator_sim(aterms)
  g <- function(nm, default) {
    if (is.null(dpars[[nm]])) rep_len(default, n) else rep_len(dpars[[nm]], n)
  }
  v <- g("mu", 0); a <- g("bs", 1); t0 <- g("ndt", 0); w <- g("bias", 0.5)
  sv <- g("sv", 0); sz <- g("sz", 0); st <- g("st", 0)
  u <- rep_len(up, n)

  nu <- numeric(n); om <- numeric(n)
  todo <- seq_len(n)
  pass <- 0L
  while (length(todo)) {
    pass <- pass + 1L
    if (pass > 1000L) {
      stop("wiener: rejection sampling for the across-trial parameters ",
           "did not converge for ", length(todo), " of ", n, " rows. ",
           "The fitted boundary probability there is essentially zero.",
           call. = FALSE)
    }
    k <- length(todo)
    nk <- stats::rnorm(k, v[todo], sv[todo])
    ok_w <- stats::runif(k, w[todo] - sz[todo] / 2, w[todo] + sz[todo] / 2)
    pu <- ddm_p_upper(nk, a[todo], ok_w)
    p <- ifelse(u[todo] == 1, pu, 1 - pu)
    keep <- stats::runif(k) < p
    nu[todo[keep]] <- nk[keep]
    om[todo[keep]] <- ok_w[keep]
    todo <- todo[!keep]
  }
  tau <- stats::runif(n, t0 - st / 2, t0 + st / 2)
  ddm_sim_rt(list(mu = nu, bs = a, ndt = tau, bias = om),
             aterms, n)
}
