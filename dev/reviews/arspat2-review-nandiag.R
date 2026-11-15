# Review of lane arspat2: where the fcor() margin route's gradient turns
# NaN (n = 400, m = 300, the failing cell of arspat2-review-nan.R). Runs
# nlminb by hand on the dry-run objective and stops at the first
# non-finite fn or gradient, printing the outer point and the pieces.
# The dry run skips mi_ml_frame(); mi_sampling(, "laplace") applies it
# and retapes, which gives the ML route frm() fits.
source("dev/reviews/arspat2-review-common.R")
k <- 20; m <- 300
W <- rook(k); W <- W / rowSums(W)
n <- nrow(W)
A0 <- diag(n) - 0.4 * W
M <- solve(crossprod(A0))
set.seed(11)
d <- data.frame(x = rnorm(n), h = factor(rep(1:8, length.out = n)))
d$y <- as.vector(solve(A0, 1 + 0.5 * d$x + rnorm(8)[d$h] +
                         rnorm(n) * sqrt(4 / rchisq(1, 4))))
set.seed(m)
d$y[sort(sample(n, m))] <- NA
pr <- set_prior("gamma(2, 0.1)", class = "nu")
ob <- quiet(frm(y | mi() ~ x + (1 | h) + fcor(M), data = d,
                data2 = list(M = M), family = student(), prior = pr,
                dry_run = "objective"))
ob <- ns$mi_sampling(ob, "laplace")
tm <- ns$autocor_seam_t_missing(ob$frame)
cat("route margin:", tm$y$margin, "m", tm$y$m, "
")
obj <- ob$obj
cat("outer:", paste(names(obj$par), collapse = " "), "\n")
trace <- list()
fn <- function(p) {
  v <- obj$fn(p)
  trace[[length(trace) + 1L]] <<- list(p = p, f = v, kind = "fn")
  v
}
gr <- function(p) {
  g <- obj$gr(p)
  trace[[length(trace) + 1L]] <<- list(p = p, g = g, kind = "gr")
  if (!all(is.finite(g))) {
    cat("non-finite gradient at p =", format(p, digits = 6), "\n")
    cat("  gradient:", format(g, digits = 4), "\n")
    rnd <- obj$env$last.par[obj$env$random]
    cat("  inner: range", format(range(rnd), digits = 4), "finite",
        all(is.finite(rnd)), "\n")
    stop("nan gradient")
  }
  g
}
t0 <- proc.time()[["elapsed"]]
r <- tryCatch(stats::nlminb(obj$par, fn, gr), error = function(e) e)
cat(sprintf("elapsed %.1f s, evaluations %d\n", proc.time()[["elapsed"]] - t0,
            length(trace)))
if (inherits(r, "error")) {
  cat("stopped:", conditionMessage(r), "\n")
  last <- utils::tail(trace, 6)
  for (t in last) {
    cat(t$kind, format(t$p, digits = 6),
        if (t$kind == "fn") format(t$f, digits = 10) else "", "\n")
  }
} else {
  cat("converged:", r$convergence, r$message, "objective", r$objective, "\n")
}
