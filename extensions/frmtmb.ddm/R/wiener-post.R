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

#' The family's `post$mean_fn`: expected response time.
#'
#' The response is the response time and the boundary is data, so the
#' mean frmtmb wants for a row is the mean response time GIVEN that
#' row's boundary: the conditional mean decision time plus the
#' non-decision time.
#'
#' @noRd
ddm_mean_rt <- function(dpars, aterms) {
  up <- aterms[["vint1"]]
  if (is.null(up)) {
    stop("wiener: the mean response time is conditional on the boundary ",
         "a trial ended at, which is missing. See ?wiener.", call. = FALSE)
  }
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
  up <- aterms[["vint1"]]
  if (is.null(up)) {
    stop("wiener: simulating a response time needs the boundary each ",
         "row ended at, which is missing. See ?wiener.", call. = FALSE)
  }
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
#' @param n Number of trials.
#' @param mu,bs,ndt,bias Drift rate, boundary separation, non-decision
#'   time and relative start point, on the response scale. Each may be a
#'   single value or a vector of length `n`, which is how a
#'   two-condition design is built.
#'
#' @return A data frame with `rt`, the response time, and `upper`, 1 for
#'   a response at the upper boundary and 0 for the lower one, in the
#'   coding `vint()` expects.
#'
#' @examples
#' set.seed(1)
#' cond <- rep(c(0, 1), each = 100)
#' dat <- ddm_simulate(200, mu = 0.2 + 1.1 * cond, bs = 1.4,
#'                     ndt = 0.3, bias = 0.5)
#' dat$cond <- factor(cond)
#' str(dat)
#'
#' @export
ddm_simulate <- function(n, mu, bs, ndt, bias = 0.5) {
  v <- rep_len(mu, n); a <- rep_len(bs, n)
  t0 <- rep_len(ndt, n); w <- rep_len(bias, n)
  if (any(a <= 0) || any(t0 < 0) || any(w <= 0) || any(w >= 1)) {
    stop("ddm_simulate(): need bs > 0, ndt >= 0 and bias strictly ",
         "inside (0, 1).", call. = FALSE)
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
