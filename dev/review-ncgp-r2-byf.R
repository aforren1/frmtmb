# Fix round 1: gp(x, by = f) against my dense marginal, with the fixed effects the fit itself holds:
# the first arm of review-ncgp-r2-struct.R gave my reference a column the model does not have.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
nug <- 1e-12
ke <- function(x, ls, lr) { K <- exp(2 * ls) * exp(-outer(x, x, "-")^2 / (2 * exp(2 * lr))); diag(K) <- diag(K) * (1 + nug); K }
for (seed in 1:3) for (sg in c(0.1, 1e-2, 1e-3)) {
  set.seed(57000 + 17 * seed + round(-log10(sg)))
  n <- 70L
  x <- sort(round(runif(n, 0, 10), 3)); z <- round(runif(n, 0, 10), 3)
  g <- factor(sample(letters[1:6], n, TRUE)); f2 <- factor(sample(c("a", "b"), n, TRUE))
  ws <- runif(n, 0.5, 2); e <- sg * rnorm(n)
  y <- (2 * sin(3 * x) + x) * (f2 == "a") + sin(x) * (f2 == "b") + e
  r <- rev_fit(y ~ gp(x, by = f2), data.frame(x, f2, y))
  X <- matrix(1, n, 1)
  a <- f2 == "a"
  f <- function(p) {
    V <- matrix(0, n, n); V[a, a] <- ke(x[a], p[2], p[3]); V[!a, !a] <- ke(x[!a], p[4], p[5]); diag(V) <- diag(V) + exp(2 * p[1])
    R <- tryCatch(chol(V), error = function(e) NULL); if (is.null(R)) return(1e300)
    wy <- backsolve(R, y, transpose = TRUE); wX <- backsolve(R, X, transpose = TRUE)
    b <- qr.solve(wX, wy); e2 <- wy - wX %*% b
    v <- 0.5 * sum(e2^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi); if (is.finite(v)) v else 1e300
  }
  sy <- sd(y); best <- Inf
  starts <- list(c(log(sigma(r$fit)[1]), r$fit$estimates$theta))
  for (ra in c(.15, .6, 2.5)) for (rb in c(.6, 2.5)) for (ns in c(.3, .01)) starts[[length(starts) + 1L]] <- c(log(sy * ns), log(sy), log(ra), log(sy), log(rb))
  for (st in starts) { o <- nlminb(st, f, control = list(iter.max = 400, rel.tol = 1e-12)); o2 <- optim(o$par, f, method = "Nelder-Mead", control = list(maxit = 500, reltol = 1e-13)); best <- min(best, o$objective, o2$value) }
  cat(sprintf("RBY s=%g seed %d names %s ll %.4f cf %.4f gap %.2e conv %d w %d\n", sg, seed, paste(names(r$fit$estimates$beta), collapse = ","), -r$fit$opt$objective, -best, -r$fit$opt$objective + best, r$fit$opt$convergence, length(r$warn)))
  flush.console()
}
cat("DONE\n")
