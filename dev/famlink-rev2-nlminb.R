## Recheck, priority 1, at the optimizer level: is nlminb on the lane's
## nlminb_trial_fn(fn) bitwise the same as nlminb on fn, when fn is NaN
## over part of the domain? No frmtmb model: synthetic objectives, so the
## mechanism is isolated. 400 random problems (seed 20260916) over:
## dimension 1 to 6; NaN region a half-space or a ball; with and without
## an analytic gradient (NaN where fn is NaN); with and without box
## bounds; start always feasible. Compared with identical() on par,
## objective, convergence, message, iterations and evaluations.
ARM <- "lane"
source("dev/famlink-rev-common.R")
wrap <- frmtmb:::nlminb_trial_fn
set.seed(20260916)
n_prob <- 400L
n_same <- 0L; n_nan <- 0L; n_warn_raw <- 0L; n_warn_wrap <- 0L
diffs <- list()
for (k in seq_len(n_prob)) {
  p <- sample(1:6, 1)
  A <- crossprod(matrix(rnorm(p * p), p)) + diag(p) * runif(1, 0.01, 2)
  target <- rnorm(p, 0, 2)
  kind <- sample(c("half", "ball", "log"), 1)
  w <- rnorm(p)
  c0 <- runif(1, -1, 1)
  bad <- switch(kind,
    half = function(x) sum(w * x) > c0 + 0.5,
    ball = function(x) sum((x - target * 0.7)^2) < 0.8,
    log = function(x) x[1] < -0.3)
  f <- function(x) {
    if (bad(x)) return(NaN)
    d <- x - target
    v <- 0.5 * sum(d * crossprod(A, d))
    if (kind == "log") v <- v - log(x[1] + 0.3)
    v
  }
  g <- function(x) {
    if (bad(x)) return(rep(NaN, p))
    d <- x - target
    gg <- as.vector(A %*% d)
    if (kind == "log") gg[1] <- gg[1] - 1 / (x[1] + 0.3)
    gg
  }
  repeat { x0 <- rnorm(p, 0, 1); if (!bad(x0) && is.finite(f(x0))) break }
  use_g <- runif(1) < 0.5
  bounded <- runif(1) < 0.4
  lo <- if (bounded) x0 - runif(p, 0.5, 3) else -Inf
  up <- if (bounded) x0 + runif(p, 0.5, 3) else Inf
  nan_seen <- 0L
  f_count <- function(x) { v <- f(x); if (is.nan(v)) nan_seen <<- nan_seen + 1L; v }
  run <- function(fn) {
    nw <- 0L
    r <- withCallingHandlers(
      tryCatch(stats::nlminb(x0, fn, if (use_g) g else NULL, lower = lo, upper = up),
               error = function(e) list(error = conditionMessage(e))),
      warning = function(wn) { nw <<- nw + 1L; invokeRestart("muffleWarning") })
    r$nwarn <- nw
    r
  }
  raw <- run(f_count)
  ww <- wrap(f)
  wr <- run(ww$fn)
  if (nan_seen > 0L) n_nan <- n_nan + 1L
  n_warn_raw <- n_warn_raw + raw$nwarn
  n_warn_wrap <- n_warn_wrap + wr$nwarn
  flds <- c("par", "objective", "convergence", "message", "iterations",
            "evaluations", "error")
  same <- vapply(flds, function(fl) identical(raw[[fl]], wr[[fl]]), TRUE)
  if (all(same)) n_same <- n_same + 1L else {
    diffs[[length(diffs) + 1L]] <- list(k = k, p = p, kind = kind, use_g = use_g,
                                        bounded = bounded, fields = flds[!same])
  }
}
cat("problems:", n_prob, " with NaN trials on the raw path:", n_nan, "\n")
cat("bitwise identical (par, objective, convergence, message, iterations, evaluations, error):",
    n_same, "of", n_prob, "\n")
cat("warnings raw:", n_warn_raw, " warnings wrapped:", n_warn_wrap, "\n")
if (length(diffs)) str(diffs[seq_len(min(10, length(diffs)))])
