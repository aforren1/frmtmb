# Review of lane nlstop, item 2569: can a cheap check detect the short
# stops of huber() with (1 | g) under the Laplace approximation? Rows of
# dev/nlstop-huber.R (100 rows, 20 groups, t3 noise, seeds 1 to 20). The
# reference: the best of the default and 8 jittered starts (as there).
# Detectors, each timed against the fit:
#   probe h: the objective at par +- h * max(1, |par|) per outer
#     coordinate (2 p inner solves); fires when one is lower by > grad_tol
#   restart: one refit from the estimate jittered by N(0, 0.3) (seed 99)
#   quad: the fit under quadrature = TRUE (the documented route)
# Usage: Rscript dev/reviews/nlstop-review-huber.R <scale>
source("dev/nlstop-common.R")
sc <- as.numeric(commandArgs(TRUE)[1])
cat("ROWHEAD scale seed ll best short silent gmax t_fit p01 t_p01 p1",
    "t_p1 rs t_rs ll_quad t_quad\n")
for (seed in 1:20) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  b <- rnorm(20)[d$g]
  d$y <- (1 + 0.5 * d$x + 0.5 * b + rt(100, 3)) * sc
  fo <- bf(y ~ x + (1 | g))
  t0 <- proc.time()[["elapsed"]]
  r <- tt(frm(fo, family = huber(), data = d))
  t_fit <- proc.time()[["elapsed"]] - t0
  if (iserr(r)) { cat("ROW", sc, seed, "ERR\n"); next }
  a <- llv(r)
  best <- a
  set.seed(seed * 31 + 7)
  for (i in 1:8) {
    st <- lapply(r$estimates[c("beta", "betad", "theta")[
      c("beta", "betad", "theta") %in% names(r$estimates)]],
      function(v) unname(v + rnorm(length(v), 0, 0.3) * pmax(1, abs(v))))
    v <- llv(tt(frm(fo, family = huber(), data = d, start = st, se = FALSE)))
    if (is.finite(v) && v > best) best <- v
  }
  short <- best - a > 1e-4
  silent <- short && !length(attr(r, "w"))
  obj <- r$obj; p0 <- r$opt$par; f0 <- r$opt$objective
  tol <- r$control$grad_tol %||% 1e-3
  gmax <- max(abs(obj$gr(p0)))
  probe <- function(h) {
    t0 <- proc.time()[["elapsed"]]
    best_gain <- -Inf
    for (i in seq_along(p0)) for (s in c(-1, 1)) {
      q <- p0; q[i] <- q[i] + s * h * max(1, abs(p0[i]))
      v <- tryCatch(obj$fn(q), error = function(e) NA_real_)
      if (is.finite(v)) best_gain <- max(best_gain, f0 - v)
    }
    obj$fn(p0)
    c(gain = best_gain, t = proc.time()[["elapsed"]] - t0)
  }
  pa <- probe(0.01); pb <- probe(0.1)
  t0 <- proc.time()[["elapsed"]]
  set.seed(99)
  st <- lapply(r$estimates[c("beta", "betad", "theta")[
    c("beta", "betad", "theta") %in% names(r$estimates)]],
    function(v) unname(v + rnorm(length(v), 0, 0.3) * pmax(1, abs(v))))
  rs <- llv(tt(frm(fo, family = huber(), data = d, start = st, se = FALSE)))
  t_rs <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  rq <- tt(frm(fo, family = huber(), data = d, quadrature = TRUE))
  t_q <- proc.time()[["elapsed"]] - t0
  cat(sprintf(paste("ROW %g %d %.6f %.6f %s %s %.3g %.2f %.3g %.2f %.3g",
                    "%.2f %.4g %.2f %.6f %.2f\n"),
              sc, seed, a, best, short, silent, gmax, t_fit,
              pa[["gain"]], pa[["t"]], pb[["gain"]], pb[["t"]], rs - a, t_rs,
              llv(rq), t_q))
  flush.console()
}
cat("DONE\n")
