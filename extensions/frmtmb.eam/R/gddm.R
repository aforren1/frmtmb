# The generalized drift-diffusion family: a Fokker-Planck solve on an RTMB tape.
#
# WHY THIS IS A ROWWISE FAMILY AND NOT A frmtmb_structure().
# The expensive object here is shared: one PDE solve serves every trial that
# shares a parameter vector, so the obvious reading is that a rowwise density
# would re-solve per row and only a structure could amortize it. That reading is
# wrong, because frmtmb calls `lpdf` ONCE per objective evaluation with
# full-length vectors, not once per row (R/objective.R, row_lpdf). The condition
# a trial belongs to is data, it arrives through vint(), and indexing an AD
# vector by data-derived integers is legal on a tape. So a rowwise density can
# do exactly one solve per condition, which is all a structure would have
# bought, and it keeps the post-fit surface: a structure that supplies loglik
# defaults every capability flag to FALSE, so fitted(), predict(), simulate()
# and residuals() would each need re-implementing or refusing.
# The two hooks that replace what frame_block() would have done are
# `aterm_data` (per-row index arithmetic, computed once from data at frame
# assembly) and `family_finalize` (constants derived from the response).

# ---------------------------------------------------------------------------
# Numerical primitives
# ---------------------------------------------------------------------------

#' Smooth positive part.
#'
#' `max(x, 0)` written as `(x + |x|)/2` because a tape cannot compare an
#' AD value against a constant, while `abs()` is overloaded.
#'
#' @noRd
gd_relu <- function(x) 0.5 * (x + abs(x))

# How many grid cells the starting mass is spread over. A delta start has
# energy at every frequency the grid carries, and Crank-Nicolson is
# A-stable but not L-stable, so the highest of those modes ring instead of
# decaying and the ringing gets worse as the spatial grid is refined at a
# fixed step. Spreading the start over two cells rather than one removes
# them. The variance it adds is of the order of the cell width squared,
# which is below the scheme's own spatial error, and it is what makes the
# density improve rather than degrade when `ny` is raised alone.
gd_start_cells <- 2

# Floors go through ddm_floor() in wiener-density.R: one helper for every
# family in the package, and the rounding trap it avoids is documented
# there.

# The floor the interpolated density is held at before it is logged.
# 1e-300 rather than the smallest positive double: this value is passed
# to log(), and log() of a denormal has a derivative near 1e320, which
# overflows to Inf and poisons the reverse pass. At 1e-300 the
# derivative stays finite, and the value is still orders of magnitude
# below anything a representable density takes, so adding it is inert
# everywhere a fit legitimately goes.
gd_dens_floor <- 1e-300

# The floor the moving boundary is held at. Chosen so that the diffusion
# coefficient it implies, 0.5 / B^2, stays many orders of magnitude below
# overflow once it is divided by the squared grid spacing; and far below
# any boundary a fitted model takes, so it is inert on every value in
# range.
gd_bound_floor <- 1e-20

#' Cubic B-spline in truncated-power form, and its integral.
#'
#' The obvious interpolation kernel, a Catmull-Rom cubic, is defined
#' piecewise on `|u|` and so branches on the very parameter it is used to
#' shift by. The cubic B-spline has a truncated-power form that needs no
#' branch, is C2 in the shift, and is a partition of unity, so a shift
#' conserves mass exactly. It smooths rather than interpolates; the
#' smoothing error is second order in the time step.
#'
#' @noRd
gd_b3 <- function(u) {
  (gd_relu(u + 2)^3 - 4 * gd_relu(u + 1)^3 + 6 * gd_relu(u)^3 -
     4 * gd_relu(u - 1)^3 + gd_relu(u - 2)^3) / 6
}

#' @noRd
gd_b3_int <- function(u) {
  (gd_relu(u + 2)^4 - 4 * gd_relu(u + 1)^4 + 6 * gd_relu(u)^4 -
     4 * gd_relu(u - 1)^4 + gd_relu(u - 2)^4) / 24
}

#' Unpack the packed tridiagonal system.
#'
#' One vector rather than four arguments because [RTMB::ADjoint()] takes a
#' function of a single vector.
#'
#' @noRd
gd_tri_unpack <- function(x) {
  n <- (length(x) + 2L) %/% 4L
  list(n = n,
       lo = x[seq_len(n - 1L)],
       di = x[(n - 1L) + seq_len(n)],
       up = x[(2L * n - 1L) + seq_len(n - 1L)],
       rhs = x[(3L * n - 2L) + seq_len(n)])
}

#' The Thomas sweep, AD-safe.
#'
#' Every accumulator is seeded from an argument so that the vectors keep
#' the AD class when the adjoint is itself taped for a second derivative.
#'
#' @noRd
gd_thomas <- function(lo, di, up, rhs) {
  "[<-" <- RTMB::ADoverload("[<-")
  n <- length(di)
  z <- di * 0
  cp <- z
  dp <- z
  cp[1L] <- up[1L] / di[1L]
  dp[1L] <- rhs[1L] / di[1L]
  for (i in 2:n) {
    m <- di[i] - lo[i - 1L] * cp[i - 1L]
    cp[i] <- if (i < n) up[i] / m else 0 * m
    dp[i] <- (rhs[i] - lo[i - 1L] * dp[i - 1L]) / m
  }
  out <- z
  out[n] <- dp[n]
  for (i in (n - 1L):1L) out[i] <- dp[i] - cp[i] * out[i + 1L]
  out
}

#' @noRd
gd_tri_f <- function(x) {
  p <- gd_tri_unpack(x)
  gd_thomas(p$lo, p$di, p$up, p$rhs)
}

#' Reverse mode of a tridiagonal solve.
#'
#' For `y = A^-1 b` and an incoming adjoint `dy`, put `w = A^-T dy`. Then
#' the adjoint of the right-hand side is `w` and the adjoint of `A` is
#' `-w_i y_j`, which is needed only on the three bands the solver stores.
#' `A^T` is the band-swapped matrix, so the transpose solve is the same
#' sweep with `lo` and `up` exchanged.
#'
#' Writing this by hand is what keeps the tape small: recorded naively the
#' sweep is a scalar recurrence of the grid's length inside every time step
#' of every condition, and RTMB records it one node at a time through an
#' R-level loop. As one atomic node the tape stops scaling with the spatial
#' resolution.
#'
#' @noRd
gd_tri_df <- function(x, y, dy) {
  p <- gd_tri_unpack(x)
  n <- p$n
  w <- gd_thomas(p$up, p$di, p$lo, dy)
  c(-w[2:n] * y[seq_len(n - 1L)], -w * y,
    -w[seq_len(n - 1L)] * y[2:n], w)
}

# The atomic is built on first use and cached: ADjoint() registers with
# RTMB, so it must happen once per session and cannot happen at build time.
gd_cache <- new.env(parent = emptyenv())

#' @noRd
gd_tri_solve <- function(x) {
  if (is.null(gd_cache$tri)) {
    gd_cache$tri <- RTMB::ADjoint(gd_tri_f, gd_tri_df, name = "gddm_tridiag")
  }
  gd_cache$tri(x)
}

#' Pick how the tridiagonal solve reaches the tape.
#'
#' Two ways, and which is faster depends on how many times the objective
#' will be evaluated, so the choice is the user's and the default is the
#' one that suits a fit.
#'
#' `"recorded"` lets RTMB record the sweep, one scalar node at a time. The
#' tape is large and slow to build, and then runs entirely in compiled
#' code. `"atomic"` collapses the sweep into a single tape node whose
#' derivative is [gd_tri_df()]. The tape is small and quick to build, but
#' every evaluation calls back into R.
#'
#' Measured at the shipped grid over six conditions, the atomic builds the
#' tape about twelve times faster and evaluates the gradient about twelve
#' times slower, so it pays only for a fit that converges in a few tens of
#' iterations. Ordinary fits do more than that, which is why `"recorded"`
#' is the default.
#'
#' @noRd
gd_tri_fn <- function(kind) {
  if (identical(kind, "atomic")) gd_tri_solve else gd_tri_f
}

#' Read one condition's value out of a per-row dpar vector.
#'
#' @noRd
gd_at <- function(v, i) if (length(v) == 1L) v[1L] else v[i]

# ---------------------------------------------------------------------------
# The solver
# ---------------------------------------------------------------------------

#' Solve the Fokker-Planck equation for one condition.
#'
#' The substitution `y = x / B(t)` pins the moving boundaries at plus and
#' minus one, so the spatial grid is fixed and no grid index depends on a
#' parameter. Under it,
#'
#'     dy = ( a(y B(t), t) / B(t) - y B'(t)/B(t) ) dt + (1 / B(t)) dW
#'
#' and every coefficient is a smooth function of every parameter. With the
#' walls stationary the scheme can be Crank-Nicolson, which is second order
#' in the step and unconditionally stable; a scheme that chases the moving
#' bound in the original coordinate cannot, because the discretized
#' operator is stiff and its eigenvalues grow as the bound collapses.
#'
#' Returns the outward probability flux at each wall on the time grid,
#' which is the defective first-passage density. The walls are absorbing,
#' so the density vanishes there and the flux is the diffusive term alone;
#' it is read with a one-sided second-order difference.
#'
#' @return A list with `up` and `lo`, each of length `nt + 1`.
#' @noRd
gd_solve <- function(p, cov, comp, ctl) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c" <- RTMB::ADoverload("c")
  nt <- ctl$nt
  ny <- ctl$ny
  dt <- ctl$dt
  h <- 2 / (ny + 1)
  yg <- seq(-1 + h, 1 - h, length.out = ny)
  i1 <- seq_len(ny - 1L)
  i2 <- 2:ny
  ones <- rep(1, ny)

  bnd <- comp$bound$fn
  drf <- comp$drift$fn
  tri <- gd_tri_fn(ctl$tridiagonal)

  rho <- comp$start$fn(yg, h, p)
  p_up <- numeric(nt + 1L) + rho[1L] * 0
  p_lo <- p_up

  # The boundary is floored for the same reason the density is, and it
  # has to be floored HERE rather than in each bound component, because
  # every component can reach zero the same way: an exponential bound
  # with a small time constant underflows to zero within the window, a
  # linear one approaches zero as its collapse fraction approaches one.
  # At B = 0 the diffusion coefficient 0.5 / B^2 is Inf, the tridiagonal
  # system is Inf on both the diagonal and the off-diagonals, and the
  # sweep returns NaN. Flooring keeps every quantity finite, so the
  # optimizer gets a bad number instead of no number in a region it only
  # ever visits by overshooting.
  b0 <- bnd(0, p, ctl)
  B0 <- ddm_floor(b0$B, gd_bound_floor)
  a0 <- drf(yg * B0, 0, p, cov) / B0 - yg * b0$dlogB
  d0 <- 0.5 / (B0 * B0)

  for (k in seq_len(nt)) {
    tk <- k * dt
    b1 <- bnd(tk, p, ctl)
    B1 <- ddm_floor(b1$B, gd_bound_floor)
    a1 <- drf(yg * B1, tk, p, cov) / B1 - yg * b1$dlogB
    d1 <- 0.5 / (B1 * B1)

    # explicit half step: rhs = (I + dt/2 M(t0)) rho, M the conservative
    # central discretization of -d/dy(a rho) + D d2/dy2 rho
    lo0 <- a0[i1] / (2 * h) + d0 / h^2
    up0 <- -a0[i2] / (2 * h) + d0 / h^2
    mv <- (-2 * d0 / h^2) * rho
    mv[i1] <- mv[i1] + up0 * rho[i2]
    mv[i2] <- mv[i2] + lo0 * rho[i1]
    rhs <- rho + (dt / 2) * mv

    # implicit half step: (I - dt/2 M(t1)) rho_new = rhs, one atomic node
    lo <- -(dt / 2) * (a1[i1] / (2 * h) + d1 / h^2)
    up <- -(dt / 2) * (-a1[i2] / (2 * h) + d1 / h^2)
    di <- ones * (1 + dt * d1 / h^2)
    new <- tri(c(lo, di, up, rhs))

    # The flux is read from the new solution alone, not averaged across
    # the step. Averaging is the natural-looking thing to do and is wrong
    # here: it values the flux at the midpoint of the step and then stores
    # it against the endpoint, which is a half-step shift of the whole
    # density. Against the analytic Wiener density that shift shows up as
    # a flat relative error of dt/2 times the density's logarithmic decay
    # rate, which does not go away as the spatial grid is refined. The new
    # solution is already second-order accurate, so reading the flux at
    # the endpoint costs no order.
    p_up[k + 1L] <- d1 * (4 * new[ny] - new[ny - 1L]) / (2 * h)
    p_lo[k + 1L] <- d1 * (4 * new[1L] - new[2L]) / (2 * h)
    rho <- new
    a0 <- a1
    d0 <- d1
  }
  list(up = p_up, lo = p_lo)
}

#' Shift a density along the time grid by the non-decision time.
#'
#' Indexing the density at `round((t - ndt)/dt)` would make an integer
#' index depend on a parameter. This is a convolution of fixed length
#' spanning the whole admissible range of `ndt`, with smooth
#' parameter-dependent weights: the indices are constants and only the
#' weights move.
#'
#' @noRd
gd_shift <- function(pv, ndt, dt, wmax) {
  "[<-" <- RTMB::ADoverload("[<-")
  n <- length(pv)
  s <- ndt / dt
  out <- numeric(n) * ndt * 0
  # The kernel is centered, so at a shift of less than two steps part of
  # it reaches backward. Those taps have to be included even though the
  # shift itself is forward: dropping them silently loses up to a sixth
  # of the mass as the non-decision time approaches zero, and the loss
  # would depend on a parameter.
  # A kernel longer than the grid has nothing left to read.
  for (j in -2L:min(wmax, n - 1L)) {
    w <- gd_b3(s - j)
    if (j >= 0L) {
      idx <- (j + 1L):n
      out[idx] <- out[idx] + w * pv[seq_len(n - j)]
    } else {
      m <- -j
      idx <- seq_len(n - m)
      out[idx] <- out[idx] + w * pv[(m + 1L):n]
    }
  }
  out
}

#' Every condition's pair of defective densities, end to end.
#'
#' The result is one long vector, condition-major and upper wall before
#' lower, so that a single gather with precomputed integer offsets reads
#' every row's density without an AD-side assignment.
#'
#' @noRd
gd_densities <- function(dpars, comp, ctl, d) {
  "c" <- RTMB::ADoverload("c")
  dt <- ctl$dt
  # Collected and concatenated once. Growing the vector with c(out, ...)
  # from an empty start does not work: the AD c() has no NULL to coerce.
  out <- vector("list", 2L * d$ncond)
  for (j in seq_len(d$ncond)) {
    i <- d$first[j]
    pj <- lapply(ctl$dpars, function(nm) gd_at(dpars[[nm]], i))
    names(pj) <- ctl$dpars
    s <- gd_solve(pj, d$cov[j, ], comp, ctl)
    pu <- gd_shift(s$up, pj$ndt, dt, ctl$wmax)
    pl <- gd_shift(s$lo, pj$ndt, dt, ctl$wmax)
    if (ctl$renormalize) {
      # The discretized solve loses probability mass, and the loss depends
      # on the parameters, because a configuration that absorbs faster
      # loses less. An unrenormalized likelihood therefore pays a hidden
      # bonus for fast absorption and biases leak and bound height.
      # Dividing by the realized mass also makes the density explicitly
      # conditional on a response inside the modeled window.
      m <- (sum(pu) + sum(pl)) * dt
      pu <- pu / m
      pl <- pl / m
    }
    if (!is.null(comp$lapse)) {
      lp <- pj$lapse
      u <- 0.5 / ctl$t_max
      pu <- (1 - lp) * pu + lp * u
      pl <- (1 - lp) * pl + lp * u
    }
    out[[2L * j - 1L]] <- pu
    out[[2L * j]] <- pl
  }
  do.call(c, out)
}

# ---------------------------------------------------------------------------
# Components
# ---------------------------------------------------------------------------

#' @noRd
gd_component <- function(kind, label, dpars, fn, aterms = character(0),
                         cov = 0L, base = FALSE) {
  structure(list(kind = kind, label = label, dpars = dpars, fn = fn,
                 aterms = aterms, cov = cov, base = base),
            class = "gddm_component")
}

#' @export
print.gddm_component <- function(x, ...) {
  cat("<gddm ", x$kind, " component: ", x$label, ">\n", sep = "")
  cat("  parameters: ", paste(names(x$dpars), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Drift specifications for [gddm()]
#'
#' The drift is the mean rate of evidence accumulation, `a(x, t)`. A GDDM
#' lets it depend on the current state `x`, on time, and on a covariate.
#' Terms are summed: pass one, or a list of several.
#'
#' \describe{
#'   \item{`gddm_drift_constant()`}{`a = mu`. The textbook drift-diffusion
#'     drift. Free parameter `mu`, identity link, so it is signed.}
#'   \item{`gddm_drift_coherence()`}{`a = sign(C) mu (|C| / cmax)^alpha`,
#'     the coherence nonlinearity of Shinn et al. (2020). The covariate `C`
#'     is supplied through `vreal()` and `cmax` scales it. Free parameters
#'     `mu` (identity link) and `alpha` (log link).}
#'   \item{`gddm_drift_leak()`}{`a = -leak x`, added to whichever base term
#'     is used. Positive `leak` is leaky integration, which pulls the
#'     accumulator back toward zero; negative `leak` is unstable
#'     integration, which pushes it away. Identity link. Note the sign is
#'     the paper's `l` and the negative of PyDDM's `leak`.}
#' }
#'
#' Exactly one base term, `gddm_drift_constant()` or
#' `gddm_drift_coherence()`, must be present and must come first, because
#' it supplies `mu`, the parameter that receives the model formula.
#'
#' @section Writing your own:
#' A drift term is the value of [gddm_drift_term()]. Its `fn` is called on
#' the tape as `fn(x, t, p, cov)` where `x` is the accumulator state at the
#' grid nodes in the original coordinate, `t` is the time as a plain
#' number, `p` is a named list holding one value of each free parameter for
#' the condition being solved, and `cov` is that condition's `vreal()`
#' values. It returns the drift at each node, and must contain no
#' comparison against a parameter.
#'
#' @param cmax Scale the coherence covariate is divided by, so that the
#'   nonlinearity is anchored at a coherence a subject sees. The paper uses
#'   the largest coherence in the design.
#' @param cov Which `vreal()` value carries the coherence, counting from 1.
#'
#' @return A `gddm_component`.
#' @seealso [gddm()], [gddm_bound_constant()], [gddm_start_point()]
#' @name gddm-drift
NULL

#' @rdname gddm-drift
#' @export
gddm_drift_constant <- function() {
  gd_component(
    "drift", "constant",
    # A drift of zero is the honest start: its sign is what the data are
    # there to say, and starting with the wrong sign costs more than
    # starting in the middle.
    list(mu = list(link = "identity", init = function(y, aterms) 0)),
    function(x, t, p, cov) x * 0 + p$mu,
    base = TRUE)
}

#' @rdname gddm-drift
#' @export
gddm_drift_coherence <- function(cmax = 1, cov = 1L) {
  if (!is.numeric(cmax) || length(cmax) != 1L || !is.finite(cmax) ||
      cmax <= 0) {
    stop("gddm_drift_coherence(): `cmax` must be one positive finite ",
         "number, the coherence the nonlinearity is anchored at.",
         call. = FALSE)
  }
  cov <- as.integer(cov)
  if (length(cov) != 1L || is.na(cov) || cov < 1L) {
    stop("gddm_drift_coherence(): `cov` must be one positive integer, ",
         "saying which vreal() value carries the coherence.", call. = FALSE)
  }
  gd_component(
    "drift", "coherence",
    list(
      # Unlike the constant drift, this one cannot start at zero: at
      # mu = 0 the drift is zero whatever alpha is, so alpha has no
      # gradient and the optimizer starts on a ridge.
      mu = list(link = "identity", init = function(y, aterms) 1),
      alpha = list(link = "log", init = function(y, aterms) 1)),
    function(x, t, p, cvv) {
      cc <- cvv[[cov]]
      # d/d(alpha) of C^alpha is C^alpha log(C), which is NaN at C = 0,
      # and a coherence design normally contains a zero condition. The
      # drift there is identically zero and carries no alpha dependence
      # at all. The coherence is data, so this branch is resolved once at
      # tape-build time and never reaches the tape.
      w <- if (cc == 0) 0 * p$mu else
        sign(cc) * p$mu * (abs(cc) / cmax)^p$alpha
      x * 0 + w
    },
    aterms = paste0("vreal", cov), cov = cov, base = TRUE)
}

#' @rdname gddm-drift
#' @export
gddm_drift_leak <- function() {
  gd_component(
    "drift", "leak",
    list(leak = list(link = "identity", init = function(y, aterms) 0)),
    function(x, t, p, cov) -p$leak * x)
}

#' Build a drift term
#'
#' The extension point behind [gddm_drift_constant()] and its siblings.
#'
#' @param label One word naming the term, used when the family prints.
#' @param dpars Named list, one entry per free parameter, each a list with
#'   `link` (a link name or a link object) and `init`, a
#'   `function(y, aterms)` returning a starting value on the natural scale.
#' @param fn `function(x, t, p, cov)` returning the drift at each node.
#' @param aterms Character vector of addition-term values `fn` reads, named
#'   as they reach the density: `"vreal1"`.
#' @param base `TRUE` if the term supplies `mu` and can stand alone.
#'
#' @return A `gddm_component`.
#' @seealso [gddm-drift]
#' @export
gddm_drift_term <- function(label, dpars, fn, aterms = character(0),
                            base = FALSE) {
  gd_component("drift", label, dpars, fn, aterms, base = base)
}

#' Boundary specifications for [gddm()]
#'
#' The decision boundaries sit at plus and minus `B(t)`, so the separation
#' between them is `2 B(t)`. All three components below name the separation
#' `bs`, as [wiener()] does, and set `B(0) = bs / 2`, so a boundary
#' separation estimated here is directly comparable with one from the
#' analytic family.
#'
#' \describe{
#'   \item{`gddm_bound_constant()`}{`B(t) = bs / 2`. Free parameter `bs`,
#'     log link.}
#'   \item{`gddm_bound_exponential()`}{`B(t) = (bs / 2) exp(-t / tau)`.
#'     Free parameters `bs` and the time constant `tau`, both log link. A
#'     large `tau` is a bound that barely collapses, which is why it is
#'     where the fit starts. Note PyDDM parameterizes the same bound by a
#'     rate, the reciprocal of `tau`.}
#'   \item{`gddm_bound_linear()`}{`B(t) = (bs / 2) (1 - kappa t / t_max)`,
#'     where `kappa` is the fraction of the bound lost by the end of the
#'     modeled window. Free parameters `bs` (log link) and `kappa` (logit
#'     link). The usual spelling of a linear collapse clips the bound at
#'     zero, which is a comparison against a parameter and cannot be taped;
#'     bounding the fraction below one instead keeps the boundary strictly
#'     positive across the window by construction, with no clipping.}
#' }
#'
#' @section Writing your own:
#' A bound term is the value of [gddm_bound_term()]. Its `fn` is called on
#' the tape as `fn(t, p, ctl)` for a plain number `t`, and returns a list
#' with `B`, the boundary at that time, and `dlogB`, its logarithmic
#' derivative `B'(t) / B(t)`. The solver needs the second because the
#' change of variable that pins the walls contributes a `-y B'(t)/B(t)`
#' term to the drift; supplying it rather than differencing `B` keeps the
#' rescaled drift exact.
#'
#' @return A `gddm_component`.
#' @seealso [gddm()], [gddm-drift], [gddm_start_point()]
#' @name gddm-bound
NULL

#' @rdname gddm-bound
#' @export
gddm_bound_constant <- function() {
  gd_component(
    "bound", "constant",
    list(bs = list(link = "log", init = function(y, aterms) 1.5)),
    function(t, p, ctl) list(B = 0.5 * p$bs, dlogB = 0 * p$bs))
}

#' @rdname gddm-bound
#' @export
gddm_bound_exponential <- function() {
  gd_component(
    "bound", "exponential",
    list(bs = list(link = "log", init = function(y, aterms) 1.5),
         tau = list(link = "log", init = function(y, aterms) max(y))),
    function(t, p, ctl) list(B = 0.5 * p$bs * exp(-t / p$tau),
                             dlogB = -1 / p$tau))
}

#' @rdname gddm-bound
#' @export
gddm_bound_linear <- function() {
  gd_component(
    "bound", "linear",
    list(bs = list(link = "log", init = function(y, aterms) 1.5),
         kappa = list(link = "logit", init = function(y, aterms) 0.1)),
    function(t, p, ctl) {
      f <- 1 - p$kappa * (t / ctl$t_max)
      list(B = 0.5 * p$bs * f, dlogB = -(p$kappa / ctl$t_max) / f)
    })
}

#' Build a boundary term
#'
#' The extension point behind [gddm_bound_constant()] and its siblings.
#'
#' @param label One word naming the term, used when the family prints.
#' @param dpars Named list, one entry per free parameter, each a list with
#'   `link` and `init`; see [gddm_drift_term()].
#' @param fn `function(t, p, ctl)` returning a list with `B` and `dlogB`.
#'
#' @return A `gddm_component`.
#' @seealso [gddm-bound]
#' @export
gddm_bound_term <- function(label, dpars, fn) {
  gd_component("bound", label, dpars, fn)
}

#' Starting-point specifications for [gddm()]
#'
#' Where the accumulator starts, as a distribution over the rescaled state
#' at time zero. `bias` is the relative start point in `(0, 1)`, the
#' fraction of the boundary separation above the lower boundary, exactly as
#' in [wiener()]; 0.5 is unbiased.
#'
#' \describe{
#'   \item{`gddm_start_point()`}{All mass at `bias`. Free parameter `bias`,
#'     logit link.}
#'   \item{`gddm_start_uniform()`}{Uniform on an interval of half-width
#'     `sz` around `bias`, with `sz` measured as a fraction of the half
#'     separation. Free parameters `bias` and `sz`, both logit link. This
#'     is the tape-safe counterpart of the Ratcliff start-point
#'     variability, which the analytic family can only reach by
#'     quadrature.}
#' }
#'
#' A point start is a delta function, which no fixed grid can hold. Both
#' components spread the initial mass with the same cubic B-spline used for
#' the non-decision-time shift, so the starting distribution is smooth in
#' `bias` and no grid index depends on a parameter. Both are then
#' normalized to unit mass on the grid, which also keeps the family honest
#' when `bias` and `sz` together push mass past a boundary.
#'
#' @section Writing your own:
#' A start term is the value of [gddm_start_term()]. Its `fn` is called as
#' `fn(y, h, p)`, where `y` are the rescaled grid nodes on `(-1, 1)` and
#' `h` is the node spacing, and it returns a density at those nodes.
#'
#' @return A `gddm_component`.
#' @seealso [gddm()], [gddm-drift], [gddm-bound]
#' @name gddm-start
NULL

#' @rdname gddm-start
#' @export
gddm_start_point <- function() {
  gd_component(
    "start", "point",
    list(bias = list(link = "logit", init = function(y, aterms) 0.5)),
    function(y, h, p) {
      hw <- h * gd_start_cells
      r <- gd_b3((y - (2 * p$bias - 1)) / hw) / hw
      r / (sum(r) * h)
    })
}

#' @rdname gddm-start
#' @export
gddm_start_uniform <- function() {
  gd_component(
    "start", "uniform",
    list(bias = list(link = "logit", init = function(y, aterms) 0.5),
         sz = list(link = "logit", init = function(y, aterms) 0.1)),
    function(y, h, p) {
      y0 <- 2 * p$bias - 1
      hw <- h * gd_start_cells
      # A smooth box: the difference of two integrated B-splines is a
      # step that is C3 in its edge, so the width differentiates cleanly.
      r <- gd_b3_int((y0 + p$sz - y) / hw) - gd_b3_int((y0 - p$sz - y) / hw)
      r / (sum(r) * h)
    })
}

#' Build a starting-point term
#'
#' The extension point behind [gddm_start_point()] and its sibling.
#'
#' @param label One word naming the term, used when the family prints.
#' @param dpars Named list, one entry per free parameter, each a list with
#'   `link` and `init`; see [gddm_drift_term()].
#' @param fn `function(y, h, p)` returning the starting density at the
#'   rescaled grid nodes.
#'
#' @return A `gddm_component`.
#' @seealso [gddm-start]
#' @export
gddm_start_term <- function(label, dpars, fn) {
  gd_component("start", label, dpars, fn)
}

# ---------------------------------------------------------------------------
# Numerical controls
# ---------------------------------------------------------------------------

#' Numerical controls for [gddm()]
#'
#' The generalized drift-diffusion likelihood has no closed form, so every
#' evaluation solves a partial differential equation on a grid. These are
#' the grid and what is done with the answer.
#'
#' @section Cost:
#' One evaluation costs one solve per **condition**, not per trial, so the
#' number of distinct parameter settings in the design is what this scales
#' with. One solve is `t_max / dt` time steps of work proportional to `ny`.
#' Halving `dt` doubles the cost and, because the scheme is second order in
#' time, divides the time-discretization error by about four. Doubling `ny`
#' doubles the cost of a step.
#'
#' The tape is built once per fit and then evaluated once per iteration.
#' At the default `tridiagonal = "recorded"` the tape grows with `t_max /
#' dt`, with `ny` and with the number of conditions, which makes the build
#' the larger cost for a short fit and the evaluation the larger cost for a
#' long one. `tridiagonal = "atomic"` moves the balance the other way.
#'
#' @section Accuracy:
#' With a constant drift and boundaries that do not move the model is the
#' Wiener model, so the solver can be checked against the closed form. At
#' the default grid, across a range of drifts, separations and starting
#' points and over decision times from 0.2 s on, the two agree to better
#' than 0.01 in the log density; a coarser grid is worse and a finer one
#' better, and the package's tests pin all three.
#'
#' Short decision times are the limit. An implicit scheme spreads a little
#' probability everywhere immediately, where the true first-passage density
#' is exponentially small, so within a few steps of zero the density is far
#' larger than the truth. It is small in absolute terms and a lapse
#' component floors it, but a fit should not be asked to read it.
#'
#' @param dt Time step, in the units of the response. The default matches
#'   PyDDM's.
#' @param ny Number of interior spatial nodes between the boundaries. Odd
#'   is the natural choice, putting an unbiased start on a node.
#' @param t_max End of the modeled window, in the units of the response.
#'   `NULL`, the default, takes the smallest multiple of `dt` strictly
#'   above the largest response time. Set it explicitly to the experiment's
#'   response deadline when there is one, and whenever you will
#'   `predict()` on new data whose largest response time differs.
#' @param renormalize Divide each condition's pair of densities by their
#'   own total mass. Leave this on. The discretized solve loses mass in a
#'   way that depends on the parameters, so a likelihood that does not
#'   renormalize rewards parameter values that absorb faster and biases the
#'   leak and the boundary height. Turning it off is provided so that the
#'   size of that bias can be measured, not because it is ever the better
#'   model.
#' @param max_ndt Upper bound for the non-decision time, in the units of
#'   the response. `NULL` takes the smallest response time. See
#'   [wiener()], whose `ndt` link this shares.
#' @param tridiagonal How the tridiagonal solve inside each step reaches
#'   the tape. `"recorded"`, the default, lets 'RTMB' record the sweep:
#'   the tape is slow to build and fast to run. `"atomic"` collapses the
#'   sweep into one node with a hand-written adjoint: the tape is quick to
#'   build and slower to run, because each evaluation calls back into R.
#'   Measured at the shipped grid over six conditions, the atomic builds
#'   about twelve times faster and evaluates about twelve times slower, so
#'   it pays only when the tape build dominates: a very fine grid, many
#'   conditions, or a fit that stops after a few tens of iterations. Both
#'   compute the same derivative.
#'
#' @return A `gddm_control`.
#' @seealso [gddm()]
#' @examples
#' gddm_control(dt = 0.005, ny = 301)
#' @export
gddm_control <- function(dt = 0.01, ny = 201L, t_max = NULL,
                         renormalize = TRUE, max_ndt = NULL,
                         tridiagonal = c("recorded", "atomic")) {
  tridiagonal <- match.arg(tridiagonal)
  if (!is.numeric(dt) || length(dt) != 1L || !is.finite(dt) || dt <= 0) {
    stop("gddm_control(): `dt` must be one positive finite number.",
         call. = FALSE)
  }
  ny <- as.integer(ny)
  if (length(ny) != 1L || is.na(ny) || ny < 5L) {
    stop("gddm_control(): `ny` must be one integer of at least 5. The ",
         "flux at each wall is read with a three-point difference, so a ",
         "grid shorter than that has nothing to read.", call. = FALSE)
  }
  if (!is.null(t_max)) {
    if (!is.numeric(t_max) || length(t_max) != 1L || !is.finite(t_max) ||
        t_max <= 0) {
      stop("gddm_control(): `t_max` must be one positive finite number, ",
           "or NULL to take it from the data.", call. = FALSE)
    }
  }
  if (!is.logical(renormalize) || length(renormalize) != 1L ||
      is.na(renormalize)) {
    stop("gddm_control(): `renormalize` must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (!is.null(max_ndt)) {
    if (!is.numeric(max_ndt) || length(max_ndt) != 1L ||
        !is.finite(max_ndt) || max_ndt <= 0) {
      stop("gddm_control(): `max_ndt` must be one positive finite ",
           "number, or NULL to take it from the data.", call. = FALSE)
    }
  }
  structure(list(dt = dt, ny = ny, t_max = t_max,
                 renormalize = renormalize, max_ndt = max_ndt,
                 tridiagonal = tridiagonal),
            class = "gddm_control")
}

#' @export
print.gddm_control <- function(x, ...) {
  cat("<gddm_control: dt = ", format(x$dt), ", ny = ", x$ny,
      ", t_max = ", if (is.null(x$t_max)) "from data" else format(x$t_max),
      ", renormalize = ", x$renormalize, ">\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# The family
# ---------------------------------------------------------------------------

#' @noRd
gd_normalize_drift <- function(drift) {
  if (inherits(drift, "gddm_component")) drift <- list(drift)
  if (!is.list(drift) || !length(drift) ||
      !all(vapply(drift, inherits, logical(1), "gddm_component"))) {
    stop("gddm(): `drift` must be one drift component or a list of them, ",
         "each from gddm_drift_constant(), gddm_drift_coherence(), ",
         "gddm_drift_leak() or gddm_drift_term().", call. = FALSE)
  }
  if (!all(vapply(drift, function(z) z$kind, character(1)) == "drift")) {
    stop("gddm(): every element of `drift` must be a drift component. A ",
         "boundary or starting-point component belongs in its own ",
         "argument.", call. = FALSE)
  }
  base <- vapply(drift, function(z) isTRUE(z$base), logical(1))
  if (sum(base) != 1L || !base[[1L]]) {
    stop("gddm(): `drift` needs exactly one base term, first in the ",
         "list. gddm_drift_constant() and gddm_drift_coherence() are ",
         "base terms; they supply `mu`, the parameter the model formula ",
         "is fitted to. gddm_drift_leak() adds to a base term and cannot ",
         "stand alone.", call. = FALSE)
  }
  drift
}

#' The generalized drift-diffusion family
#'
#' The generalized drift-diffusion model of Shinn, Lam and Murray (2020).
#' Evidence accumulates between two decision boundaries and the response
#' time is the first time it reaches one, as in [wiener()]; unlike
#' [wiener()], the drift may depend on the state of the accumulator and on
#' a covariate, and the boundaries may move. That generality costs the
#' closed form: there is no first-passage density to evaluate, so each
#' likelihood evaluation solves the Fokker-Planck equation forward in time
#' and reads the probability flux through each boundary.
#'
#' @section The model:
#' The accumulator `x` follows
#'
#' \deqn{dx = a(x, t) \, dt + dW}
#'
#' from a starting distribution `X_0`, absorbed at moving boundaries
#' \eqn{\pm B(t)}. The response time is the absorption time plus a
#' non-decision time `ndt`. The drift `a`, the boundary `B` and the start
#' `X_0` are chosen by argument; each brings its own free parameters, and
#' each of those takes a formula like any other distributional parameter.
#' See [gddm-drift], [gddm-bound] and [gddm-start] for the catalogue and
#' for how to add to it.
#'
#' The likelihood is solved after the substitution \eqn{y = x / B(t)},
#' which pins the boundaries at \eqn{\pm 1} and leaves the grid fixed while
#' the boundary collapses. That is what makes the model differentiable: the
#' usual treatment sandwiches a moving bound between two integer grid
#' indices, which makes the objective a function of where the bound falls
#' between nodes and rules out gradient-based fitting. With the walls
#' stationary the scheme is Crank-Nicolson.
#'
#' @section What the data must carry:
#' Which boundary a trial ended at is data, and so is the condition a trial
#' belongs to. Two spellings carry the boundary, as they do for
#' [wiener()]:
#'
#' ```
#' frm(bf(rt | dec(response) + vint(cond) ~ 1), family = gddm(), data = dat)
#' frm(bf(rt | vint(upper, cond) ~ 1),          family = gddm(), data = dat)
#' ```
#'
#' `dec()` takes what brms takes, a factor whose second level is the upper
#' boundary or a 0/1 column. `upper` in the `vint()` spelling is the same
#' 0/1. Two is all there is: one accumulator in one dimension between two
#' absorbing boundaries can end a trial at the upper wall or the lower
#' wall and nowhere else, so multi-alternative choice is outside this
#' family, and a third level is refused and pointed at `lba()` rather than
#' folded into one of the two.
#'
#' `cond` is an integer labelling the distinct parameter settings in the
#' design: one solve serves every trial that shares one. Note where it
#' sits. `vint()` numbers its values positionally, so the condition is the
#' first `vint()` value when `dec()` carries the boundary and the second
#' when `vint()` carries it. The family reads whichever it is; you write
#' the pair in the order the two lines above show. A drift term that reads
#' a covariate needs `vreal()` as well, for example
#' `bf(rt | dec(response) + vint(cond) + vreal(coh) ~ 1)`.
#'
#' **Every row sharing a condition must share every parameter value.** The
#' family cannot check this, because checking would mean comparing values
#' on the tape, so it is your side of the contract: build the index from
#' every variable that appears on the right-hand side of any of the
#' family's formulas, which is what [gddm_conditions()] does. What the
#' family can check, and does, is that the `vreal()` covariates are
#' constant within a condition.
#'
#' @section Cost, honestly:
#' This is much slower than [wiener()] and you should use [wiener()]
#' whenever it applies, which is whenever the drift is constant in the
#' state and in time and the boundaries do not move. There the
#' first-passage density is a known pair of series and costs arithmetic;
#' here every evaluation solves a PDE per condition. What you buy is the
#' class of models the analytic density cannot express at all: collapsing
#' boundaries, leaky or unstable integration, a drift that varies within a
#' trial. Cost scales with the number of conditions, so a design with many
#' distinct parameter settings is where this becomes painful.
#' See [gddm_control()].
#'
#' @param drift A drift component, or a list of them to be summed. The
#'   first must be a base term. See [gddm-drift].
#' @param bound A boundary component. See [gddm-bound].
#' @param start A starting-point component. See [gddm-start].
#' @param lapse `"none"`, or `"uniform"` to mix the first-passage density
#'   with a lapse distribution uniform over the modeled window and split
#'   evenly between the two responses. `"uniform"` adds the free parameter
#'   `lapse` on a logit link; the published models fix it rather than
#'   fitting it, which is spelled `bf(..., lapse = 0.05)`.
#' @param control Numerical controls, from [gddm_control()].
#'
#' @return A `frmtmb_family`.
#'
#' @references
#' Shinn, M., Lam, N. H. and Murray, J. D. (2020). A flexible framework for
#' simulating and fitting generalized drift-diffusion models. *eLife*, 9,
#' e56938.
#'
#' @seealso [wiener()] for the analytic special case, [gddm_control()] for
#'   the grid, [gddm_conditions()] for building the condition index, and
#'   `vignette("ddm")`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' # a coarse grid, so that the example is quick; see gddm_control()
#' ctl <- gddm_control(t_max = 2, dt = 0.02, ny = 101)
#' dat <- gddm_simulate(400, mu = 2.5, bs = 3, ndt = 0.25, tau = 1,
#'                      bound = gddm_bound_exponential(), control = ctl)
#' dat$cond <- 1L
#' fit <- frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
#'            family = gddm(bound = gddm_bound_exponential(), control = ctl),
#'            data = dat)
#' fixef(fit)
#' }
#' @export
gddm <- function(drift = gddm_drift_constant(),
                 bound = gddm_bound_constant(),
                 start = gddm_start_point(),
                 lapse = c("none", "uniform"),
                 control = gddm_control()) {
  drift <- gd_normalize_drift(drift)
  if (!inherits(bound, "gddm_component") || bound$kind != "bound") {
    stop("gddm(): `bound` must be one boundary component, from ",
         "gddm_bound_constant(), gddm_bound_exponential(), ",
         "gddm_bound_linear() or gddm_bound_term().", call. = FALSE)
  }
  if (!inherits(start, "gddm_component") || start$kind != "start") {
    stop("gddm(): `start` must be one starting-point component, from ",
         "gddm_start_point(), gddm_start_uniform() or ",
         "gddm_start_term().", call. = FALSE)
  }
  if (!inherits(control, "gddm_control")) {
    stop("gddm(): `control` must come from gddm_control().", call. = FALSE)
  }
  lapse <- match.arg(lapse)

  # One summed drift closure, built once, so the tape sees a single
  # expression rather than a dispatch per node.
  dfns <- lapply(drift, function(z) z$fn)
  drift_fn <- function(x, t, p, cov) {
    v <- dfns[[1L]](x, t, p, cov)
    for (k in seq_along(dfns)[-1L]) v <- v + dfns[[k]](x, t, p, cov)
    v
  }
  comp <- list(
    drift = list(fn = drift_fn,
                 label = paste(vapply(drift, function(z) z$label,
                                      character(1)), collapse = " + "),
                 aterms = unique(unlist(lapply(drift,
                                               function(z) z$aterms)))),
    bound = bound, start = start,
    lapse = if (lapse == "none") NULL else lapse)

  # The non-decision time is always free and always last, so that the
  # solver's parameters come first and read in the order of the model.
  terms <- c(drift, list(bound, start))
  dp <- list()
  for (tm in terms) {
    for (nm in names(tm$dpars)) {
      if (!is.null(dp[[nm]])) {
        stop("gddm(): two components both supply the parameter `", nm,
             "`. Each free quantity must be named once.", call. = FALSE)
      }
      dp[[nm]] <- tm$dpars[[nm]]
    }
  }
  # The real link is a logit scaled onto (0, max_ndt), and the bound comes
  # from the data. Until family_finalize() supplies it the family says so
  # rather than naming a link it will not use.
  pending <- function(...) {
    stop("gddm: the non-decision-time link is not resolved yet. Its ",
         "upper bound comes from the response, so it is set when frm() ",
         "assembles the model frame. Pass max_ndt to gddm_control() to ",
         "fix the bound up front and inspect the family before a fit.",
         call. = FALSE)
  }
  dp[["ndt"]] <- list(
    link = list(name = "scaled_logit(0, from data)", linkfun = pending,
                linkinv = pending, mu_eta = pending),
    init = function(y, aterms) 0.5 * min(y))
  if (lapse == "uniform") {
    dp[["lapse"]] <- list(link = "logit", init = function(y, aterms) 0.02)
  }
  if (names(dp)[[1L]] != "mu") {
    stop("gddm(): the first free parameter must be `mu`. A base drift ",
         "term supplies it; a component that renames it cannot receive ",
         "the model formula.", call. = FALSE)
  }
  dpnames <- names(dp)
  # Only the covariates can be declared. Neither the boundary nor the
  # condition can: the boundary arrives as dec() OR vint1, and the
  # condition's own slot moves with that choice, and a declaration
  # cannot say "either". Both are refused by hand in gd_check_response().
  req <- unique(unlist(lapply(terms, function(z) z$aterms)))
  if (is.null(req)) req <- character(0)

  fam <- frmtmb::custom_family(
    "gddm",
    dpars = dpnames,
    links = lapply(dp, function(z) z$link),
    lpdf = function(y, dpars, aterms) {
      stop("gddm: this family was used without being finalized against ",
           "the data. That happens only if the family object is called ",
           "outside frm(); the grid and the non-decision-time bound are ",
           "resolved when the model frame is assembled.", call. = FALSE)
    },
    valid_y = function(y, aterms) gd_check_response(y, aterms, comp),
    init_dpars = lapply(dp, function(z) z$init),
    type = "continuous",
    required_aterms = req,
    family_finalize = function(fam, y, aterms) {
      gd_finalize(fam, y, aterms, comp, control, dpnames)
    },
    sim_refusal = NULL)
  fam[["gddm"]] <- list(comp = comp, control = control, dpars = dpnames)
  fam
}

#' Response, boundary and condition validation.
#'
#' The decision indicator and the condition index, under whichever
#' spelling supplied them.
#'
#' `dec()` is the addition term this package contributes to frmtmb's
#' registry and the one every reference on the model uses; `vint()` is
#' the general-purpose route that was the only one available before that
#' registry existed, and it keeps working. `wiener()` reads its boundary
#' the same way.
#'
#' The condition index moves with the spelling, and it has to. `vint()`
#' numbers its values positionally, so when `dec()` carries the boundary
#' the user's `vint(cond)` is `vint1`, not `vint2`. Reading a fixed slot
#' would make a `dec()` user write a dummy first value to push the
#' condition into place.
#'
#' @noRd
gd_indicator <- function(aterms) {
  d <- aterms[["dec"]]
  if (is.null(d)) {
    list(up = aterms[["vint1"]], cond = aterms[["vint2"]],
         cond_is = "the second value of vint()")
  } else {
    list(up = d, cond = aterms[["vint1"]],
         cond_is = "the first value of vint()")
  }
}

#' Runs once when the model frame is assembled, before any link is used.
#'
#' @noRd
gd_check_response <- function(y, aterms, comp) {
  if (any(!is.finite(y)) || any(y <= 0)) {
    stop("gddm: the response must be a strictly positive, finite ",
         "response time.", call. = FALSE)
  }
  ix <- gd_indicator(aterms)
  up <- ix[["up"]]
  if (is.null(up)) {
    # `required_aterms` names the terms a density needs ALL of, and this
    # family needs EITHER of two spellings for the boundary, so the
    # refusal is written out here instead. wiener() carries the same
    # hand-rolled check for the same reason.
    stop("gddm: the decision indicator is missing. Which boundary a ",
         "trial ended at is data, and it reaches the family through ",
         "dec(), as it does in brms:
",
         "    frm(bf(rt | dec(response) + vint(cond) ~ x), ",
         "family = gddm(), ...)
",
         "where `response` is a factor whose second level is the upper ",
         "boundary, or a 0/1 column. vint(upper, cond) carries the ",
         "same pair as plain integers and also works, boundary first ",
         "and condition second.", call. = FALSE)
  }
  if (any(!is.finite(up))) {
    stop("gddm: the decision indicator holds a missing or infinite ",
         "value. Which boundary a trial ended at is data and has to be ",
         "known for every trial the density scores.", call. = FALSE)
  }
  # The model is one accumulator between two absorbing boundaries, so the
  # number of responses it can express is two, structurally. A third
  # level is not a parameterization of this model, and silently folding
  # it into one of the two boundaries would fit something nobody asked
  # for, so it is refused here, at the first point the family sees data.
  lev <- sort(unique(up))
  if (length(lev) > 2L) {
    stop("gddm: the decision indicator has ", length(lev), " distinct ",
         "values, and this family admits exactly two. A generalized ",
         "drift-diffusion model is a single accumulator between two ",
         "absorbing boundaries, so a trial can end at the upper one or ",
         "the lower one and nowhere else. More than two alternatives is ",
         "a different architecture, not another parameter: it needs ",
         "racing accumulators rather than one accumulator between two ",
         "walls. That is what lba() in this package fits, so use ",
         "lba(n) for n alternatives. Collapse the response to two if ",
         "you want this family.",
         call. = FALSE)
  }
  if (!all(lev %in% c(0, 1))) {
    stop("gddm: the decision indicator must be 0 at the lower boundary ",
         "and 1 at the upper one. dec() reads a factor on its levels ",
         "and produces that coding for you, taking the SECOND level as ",
         "the upper boundary; a numeric column is passed through as it ",
         "stands, so recode one yourself with ",
         "as.integer(decision == \"upper\").", call. = FALSE)
  }
  cnd <- ix[["cond"]]
  if (is.null(cnd)) {
    stop("gddm: the condition index is missing. One solve of the ",
         "Fokker-Planck equation serves every trial that shares a ",
         "parameter vector, and this family finds those trials through ",
         "an index it is given, because comparing parameter values is ",
         "not something a tape can do. Supply it as ", ix[["cond_is"]],
         ", which gddm_conditions() builds.", call. = FALSE)
  }
  if (any(!is.finite(cnd)) || any(cnd < 1)) {
    stop("gddm: the condition index must be a positive integer ",
         "labelling the distinct parameter settings in the design. It ",
         "is ", ix[["cond_is"]], " here, and gddm_conditions() builds ",
         "one.", call. = FALSE)
  }
  # The covariates a drift term reads are the one part of the
  # constant-within-condition contract the family can check, so it does.
  for (nm in unique(unlist(lapply(comp, function(z)
    if (is.list(z)) z$aterms else NULL)))) {
    v <- aterms[[nm]]
    if (is.null(v)) next
    ok <- tapply(v, cnd, function(z) length(unique(z)) == 1L)
    if (!all(ok)) {
      stop("gddm: ", nm, " is not constant within every condition. One ",
           "solve serves a whole condition, so a covariate the drift ",
           "reads has to take one value there. Split the condition ",
           "index, or build it with gddm_conditions().", call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Resolve the grid, the links and the density against the data.
#'
#' `family_finalize` is the sanctioned seam for a family that derives
#' something from the response: whatever it returns is the family the rest
#' of the fit sees. Everything downstream of it, the density and the
#' per-row index arithmetic alike, is rebuilt here rather than reading a
#' mutable environment, so the family object never lies about what it is.
#'
#' @noRd
gd_finalize <- function(fam, y, aterms, comp, control, dpnames) {
  # Idempotent on purpose. predict() and simulate() on new data reassemble
  # the frame, and a grid re-derived from new response times would be a
  # different model from the one that was fitted.
  if (!is.null(fam[["gddm"]][["ctl"]])) return(fam)
  dt <- control$dt
  t_max <- control$t_max
  if (is.null(t_max)) {
    t_max <- (floor(max(y) / dt) + 1) * dt
  } else if (t_max <= max(y)) {
    stop("gddm: t_max = ", format(t_max), " is at or below the largest ",
         "response time (", format(max(y)), "). The modeled window has ",
         "to contain every response it is asked to score.", call. = FALSE)
  }
  ub <- control$max_ndt
  if (is.null(ub)) {
    ub <- min(y)
  } else if (ub > min(y)) {
    stop("gddm: max_ndt = ", format(ub), " is above the smallest ",
         "response time (", format(min(y)), "). The density is zero at ",
         "and below the non-decision time, so a bound above min(rt) ",
         "admits parameter values with no likelihood.", call. = FALSE)
  }
  nt <- as.integer(round(t_max / dt))
  ctl <- list(dt = dt, ny = control$ny, t_max = nt * dt, nt = nt,
              renormalize = control$renormalize, max_ndt = ub,
              wmax = as.integer(ceiling(ub / dt)) + 2L, dpars = dpnames,
              tridiagonal = control$tridiagonal)

  # The same scaled logit wiener() uses: the density is zero at and below
  # the non-decision time, so the constraint is made structural rather
  # than left for the optimizer to discover by walking over a cliff.
  fam[["links"]][["ndt"]] <- list(
    name = paste0("scaled_logit(0, ", signif(ub, 4), ")"),
    linkfun = function(mu) log(mu / (ub - mu)),
    linkinv = function(eta) ub / (1 + exp(-eta)),
    mu_eta = function(eta) {
      p <- 1 / (1 + exp(-eta))
      ub * p * (1 - p)
    })

  fam[["aterm_data"]] <- function(y, aterms) {
    gd_aterm_data(y, aterms, comp, ctl)
  }
  # A response time is data, so which two grid nodes it falls between and
  # the weight between them are constants: the density is one gather.
  fam[["lpdf"]] <- frmtmb::frmtmb_ad_overload(function(y, dpars, aterms) {
    d <- aterms[[".gddm"]]
    p <- gd_densities(dpars, comp, ctl, d)
    # Floored, not left to underflow. The optimizer walks through
    # parameter values whose density at some observed response time the
    # grid cannot represent, and an unfloored log() answers NaN there.
    # NaN is not a value an optimizer can use, and inside mixture() it
    # takes every other component with it through the log-sum-exp. See
    # gddm_floored() for how many rows this caught at the optimum.
    log(ddm_floor((1 - d$w) * p[d$i1] + d$w * p[d$i2], gd_dens_floor))
  })
  fam[["post"]] <- list(mean_fn = function(dpars, aterms) {
    gd_mean_rt(dpars, aterms, comp, ctl)
  })
  fam[["sim"]] <- function(dpars, aterms, n) {
    gd_sim_rt(dpars, aterms, comp, ctl, n)
  }
  fam[["gddm"]] <- list(comp = comp, control = control, ctl = ctl,
                        dpars = dpnames)
  fam
}

#' Per-row index arithmetic, done once from data.
#'
#' Response times are data, so the two grid nodes a response falls between
#' and the weight between them are constants. Precomputing them here means
#' the density does one gather from one long vector instead of an AD-side
#' assignment per condition, and it is the rowwise family's answer to what
#' a structure would have used `frame_block()` for.
#'
#' @noRd
gd_aterm_data <- function(y, aterms, comp, ctl) {
  ix <- gd_indicator(aterms)
  cnd <- as.integer(ix[["cond"]])
  lev <- sort(unique(cnd))
  j <- match(cnd, lev)
  ncond <- length(lev)
  first <- match(seq_len(ncond), j)

  ncov <- 0L
  for (z in comp) {
    if (is.list(z) && length(z$aterms)) {
      ncov <- max(ncov, max(as.integer(sub("^vreal", "", z$aterms))))
    }
  }
  cov <- matrix(0, nrow = ncond, ncol = max(ncov, 1L))
  if (ncov > 0L) {
    for (k in seq_len(ncov)) {
      v <- aterms[[paste0("vreal", k)]]
      if (!is.null(v)) cov[, k] <- v[first]
    }
  }

  s <- y / ctl$dt
  k0 <- floor(s)
  w <- s - k0
  k0 <- as.integer(k0)
  nb <- ctl$nt + 1L
  # The window is fixed when the model is fitted, so new data can carry a
  # response time the grid was never built to reach. Caught here rather
  # than read off the end of the density.
  if (any(k0 < 0L) || any(k0 + 2L > nb)) {
    stop("gddm: a response time of ", format(max(y)), " falls outside the ",
         "modeled window, which ends at ", format(ctl$t_max), ". The ",
         "window is fixed when the model is fitted, so set t_max in ",
         "gddm_control() high enough to cover the data you will predict ",
         "on as well as the data you fit to.", call. = FALSE)
  }
  up <- as.integer(ix[["up"]])
  # condition-major, upper wall then lower, so one offset reaches any row
  base <- (j - 1L) * 2L * nb + (1L - up) * nb
  list(.gddm = list(ncond = ncond, first = first, cov = cov, gindex = j,
                    w = w, i1 = base + k0 + 1L, i2 = base + k0 + 2L))
}

#' Draw times from a density given on the time grid.
#'
#' The distribution function is built by the trapezoid rule at the nodes,
#' which is the distribution function of the piecewise-linear density the
#' likelihood reads when it interpolates between two nodes. Summing the
#' nodes instead would put the mass of each cell half a step late, and the
#' draws would then come from a slightly different model than the one being
#' fitted: measured over replicates, that offset biased the drift rate by
#' several Monte Carlo standard errors while leaving the other parameters
#' alone.
#'
#' @noRd
gd_draw <- function(v, tg, n) {
  v[!is.finite(v) | v < 0] <- 0
  m <- length(v)
  cdf <- c(0, cumsum((v[-1L] + v[-m]) / 2))
  if (cdf[m] <= 0) return(rep(NA_real_, n))
  stats::approx(cdf / cdf[m], tg, xout = stats::runif(n),
                ties = "ordered", rule = 2)$y
}

#' Conditional mean response time, for fitted() and predict().
#'
#' Runs at the estimates, off the tape, so plain numbers throughout.
#'
#' @noRd
gd_mean_rt <- function(dpars, aterms, comp, ctl) {
  d <- aterms[[".gddm"]]
  p <- as.numeric(gd_densities(dpars, comp, ctl, d))
  tg <- seq(0, ctl$t_max, by = ctl$dt)
  nb <- ctl$nt + 1L
  wall <- gd_row_wall(d, nb)
  out <- rep(NA_real_, length(d$gindex))
  for (j in seq_len(d$ncond)) {
    for (b in 0:1) {
      rows <- which(d$gindex == j & wall == b)
      if (!length(rows)) next
      v <- p[(j - 1L) * 2L * nb + b * nb + seq_len(nb)]
      m <- sum(v)
      if (m > 0) out[rows] <- sum(tg * v) / m
    }
  }
  out
}

#' Which wall each row read, recovered from its own gather offset.
#'
#' 0 is the upper boundary and 1 the lower, the order the densities are
#' concatenated in.
#'
#' @noRd
gd_row_wall <- function(d, nb) ((d$i1 - 1L) %/% nb) %% 2L

#' Draw a response time for each row, conditional on its own boundary.
#'
#' @noRd
gd_sim_rt <- function(dpars, aterms, comp, ctl, n) {
  d <- aterms[[".gddm"]]
  p <- as.numeric(gd_densities(dpars, comp, ctl, d))
  tg <- seq(0, ctl$t_max, by = ctl$dt)
  nb <- ctl$nt + 1L
  wall <- gd_row_wall(d, nb)
  out <- numeric(n)
  for (j in seq_len(d$ncond)) {
    for (b in 0:1) {
      rows <- which(d$gindex == j & wall == b)
      if (!length(rows)) next
      out[rows] <- gd_draw(p[(j - 1L) * 2L * nb + b * nb + seq_len(nb)],
                           tg, length(rows))
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Helpers a user needs
# ---------------------------------------------------------------------------

#' Build a condition index for [gddm()]
#'
#' One solve of the generalized drift-diffusion likelihood serves every
#' trial that shares a parameter vector, and the family finds those trials
#' through an index it is given rather than by comparing parameter values,
#' which it cannot do on a tape. This builds the index: the distinct
#' combinations of the variables you name, numbered.
#'
#' Name every variable that appears on the right-hand side of any formula
#' in the model, and every covariate a drift term reads. Naming more than
#' that is safe and only costs solves; naming fewer is wrong, and wrong in
#' a way nothing downstream can detect.
#'
#' @param data A data frame.
#' @param ... Bare variable names, or a one-sided formula naming them.
#'
#' @return An integer vector, one entry per row of `data`.
#' @seealso [gddm()]
#' @examples
#' d <- data.frame(coh = c(0, 0, 0.5, 0.5), block = c(1, 2, 1, 2))
#' gddm_conditions(d, coh, block)
#' gddm_conditions(d, ~ coh)
#' @export
gddm_conditions <- function(data, ...) {
  if (!is.data.frame(data)) {
    stop("gddm_conditions(): `data` must be a data frame.", call. = FALSE)
  }
  # The dots are evaluated in `data` first: a bare variable name is a
  # column, and evaluating it in the caller's frame would not find one.
  exprs <- as.list(substitute(list(...)))[-1L]
  vals <- lapply(exprs, eval, data, parent.frame())
  if (length(vals) == 1L && inherits(vals[[1L]], "formula")) {
    vars <- all.vars(vals[[1L]])
    miss <- setdiff(vars, names(data))
    if (length(miss)) {
      stop("gddm_conditions(): the formula names ",
           paste(miss, collapse = ", "), ", which `data` does not have.",
           call. = FALSE)
    }
    cols <- data[vars]
  } else {
    cols <- as.data.frame(vals)
  }
  if (!ncol(cols)) {
    stop("gddm_conditions(): name at least one variable. An index with ",
         "no variables in it puts every trial in one condition, which is ",
         "right only when nothing in the model varies across rows.",
         call. = FALSE)
  }
  as.integer(factor(do.call(paste, c(cols, sep = "\r"))))
}

#' How many rows the grid could not represent
#'
#' The generalized likelihood is floored: where the solved density at a
#' trial's own response time underflows, the log density is a large
#' finite negative number rather than `NaN`, so the optimizer gets a
#' value it can use and a mixture's log-sum-exp is not poisoned by one
#' component. That is the right behavior and it is silent, which is why
#' this exists: it says, once and after the fit, how many rows were
#' answered by the floor rather than by the solver.
#'
#' A count of zero is the ordinary case and means the grid represented
#' every observation. A small count means a few trials sit in the
#' leading edge, within a few time steps of the fitted non-decision
#' time, where a fixed grid cannot resolve a first-passage density that
#' is climbing through orders of magnitude; those rows contributed a
#' constant instead of information. A large count means the fit is not
#' to be trusted: shrink `dt` in [gddm_control()], or give the model a
#' lapse component with `gddm(lapse = "uniform")`, which floors the
#' density in the model rather than in the arithmetic.
#'
#' @param fit A fitted `gddm()` model.
#'
#' @return The number of rows at the floor, with the row indices in the
#'   `"rows"` attribute and the number of observations in `"n_obs"`.
#' @seealso [gddm()], [gddm_control()]
#' @examples
#' \donttest{
#' set.seed(3)
#' ctl <- gddm_control(t_max = 2, dt = 0.02, ny = 101)
#' dat <- gddm_simulate(200, mu = 2, bs = 2.5, ndt = 0.25, control = ctl)
#' dat$cond <- 1L
#' fit <- frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
#'            family = gddm(control = ctl), data = dat)
#' gddm_floored(fit)
#' }
#' @export
gddm_floored <- function(fit) {
  rsp <- frmtmb::single_response(fit)
  famobj <- rsp[["family"]]
  if (!identical(famobj[["family"]], "gddm")) {
    stop("gddm_floored(): this is a ", famobj[["family"]], " model. The ",
         "floor it reports belongs to the generalized drift-diffusion ",
         "density, so there is nothing to count on a fit from another ",
         "family.", call. = FALSE)
  }
  bag <- famobj[["gddm"]]
  atv <- fit$frame[["aterm_values"]][[rsp[["resp_name"]]]]
  gd <- atv[[".gddm"]]
  dvals <- frmtmb::eval_dpars(fit)[[rsp[["resp_name"]]]]
  pv <- as.numeric(gd_densities(dvals, bag[["comp"]], bag[["ctl"]], gd))
  # the raw interpolated density, before the floor the likelihood adds
  raw <- (1 - gd$w) * pv[gd$i1] + gd$w * pv[gd$i2]
  rows <- which(!(raw > gd_dens_floor))
  structure(length(rows), rows = rows, n_obs = length(raw))
}

#' Simulate from a generalized drift-diffusion model
#'
#' Draws choices and response times by solving the model's own
#' Fokker-Planck equation and sampling the resulting defective densities,
#' so the draws come from the density [gddm()] fits rather than from a
#' discretized forward simulation, whose first passages are late by however
#' much the step misses excursions between monitoring times.
#'
#' @param n Number of trials.
#' @param ... Parameter values by name: `mu`, `bs`, `ndt` and whatever else
#'   the chosen components need (`alpha`, `leak`, `tau`, `kappa`, `bias`,
#'   `sz`, `lapse`). Anything not given takes the component's own starting
#'   value.
#' @param coh Coherence covariate, recycled to length `n`. Only meaningful
#'   with [gddm_drift_coherence()].
#' @param drift,bound,start,lapse,control As in [gddm()].
#'
#' @return A data frame with `rt`, `upper`, `cond` and, when a coherence
#'   drift is used, `coh`.
#' @seealso [gddm()]
#' @examples
#' set.seed(2)
#' head(gddm_simulate(20, mu = 1.5, bs = 2, ndt = 0.2,
#'                    control = gddm_control(t_max = 2, dt = 0.02,
#'                                           ny = 101)))
#' @export
gddm_simulate <- function(n, ..., coh = 0,
                          drift = gddm_drift_constant(),
                          bound = gddm_bound_constant(),
                          start = gddm_start_point(),
                          lapse = c("none", "uniform"),
                          control = gddm_control()) {
  lapse <- match.arg(lapse)
  fam <- gddm(drift = drift, bound = bound, start = start, lapse = lapse,
              control = control)
  comp <- fam[["gddm"]]$comp
  dpnames <- fam[["gddm"]]$dpars
  vals <- list(...)
  bad <- setdiff(names(vals), dpnames)
  if (length(bad)) {
    stop("gddm_simulate(): the chosen components have no parameter ",
         paste(bad, collapse = ", "), ". They take ",
         paste(dpnames, collapse = ", "), ".", call. = FALSE)
  }
  coh <- rep_len(coh, n)
  cnd <- as.integer(factor(coh))
  ncond <- max(cnd)
  first <- match(seq_len(ncond), cnd)

  dt <- control$dt
  t_max <- if (is.null(control$t_max)) 2 else control$t_max
  nt <- as.integer(round(t_max / dt))

  defaults <- list(mu = 1, alpha = 1, leak = 0, bs = 1.5, tau = t_max,
                   kappa = 0.1, bias = 0.5, sz = 0.1, ndt = 0.2,
                   lapse = 0.02)
  pv <- lapply(dpnames, function(nm) {
    v <- vals[[nm]]
    if (is.null(v)) v <- defaults[[nm]]
    if (is.null(v)) {
      stop("gddm_simulate(): no value and no default for `", nm,
           "`. Give it by name.", call. = FALSE)
    }
    rep_len(as.numeric(v), n)
  })
  names(pv) <- dpnames
  if (max(pv$ndt) >= t_max) {
    stop("gddm_simulate(): ndt is at or past the end of the simulated ",
         "window, so no trial can produce a response time inside it. ",
         "Raise t_max in control, or lower ndt.", call. = FALSE)
  }
  # The shift kernel only has to span the non-decision times actually
  # asked for, unlike a fit, where it has to span everything the link
  # can reach.
  ctl <- list(dt = dt, ny = control$ny, t_max = nt * dt, nt = nt,
              renormalize = control$renormalize, max_ndt = max(pv$ndt),
              wmax = as.integer(ceiling(max(pv$ndt) / dt)) + 2L,
              dpars = dpnames, tridiagonal = control$tridiagonal)
  d <- list(ncond = ncond, first = first,
            cov = matrix(coh[first], ncol = 1L), gindex = cnd)
  p <- gd_densities(pv, comp, ctl, d)
  tg <- seq(0, ctl$t_max, by = dt)
  nb <- nt + 1L
  rt <- numeric(n)
  upper <- integer(n)
  for (j in seq_len(ncond)) {
    rows <- which(cnd == j)
    vu <- as.numeric(p[(j - 1L) * 2L * nb + seq_len(nb)])
    vl <- as.numeric(p[(j - 1L) * 2L * nb + nb + seq_len(nb)])
    vu[!is.finite(vu) | vu < 0] <- 0
    vl[!is.finite(vl) | vl < 0] <- 0
    pu <- sum(vu) / (sum(vu) + sum(vl))
    u <- stats::runif(length(rows)) < pu
    upper[rows] <- as.integer(u)
    for (b in c(TRUE, FALSE)) {
      sel <- rows[u == b]
      if (!length(sel)) next
      rt[sel] <- gd_draw(if (b) vu else vl, tg, length(sel))
    }
  }
  out <- data.frame(rt = rt, upper = upper, cond = cnd)
  if (length(unlist(lapply(comp, function(z)
    if (is.list(z)) z$aterms else NULL)))) {
    out$coh <- coh
  }
  out
}
