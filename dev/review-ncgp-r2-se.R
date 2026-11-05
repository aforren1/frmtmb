# Fix round 1, item 6, 11: the five fits where se_row_curved() changes a gp fit at ordinary noise, read
# against the closed form's Hessian. Is the rescued SE right? Which parameter? Seeds as review-ncgp-r2-claim.R.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
patch <- function(nm, f) { environment(f) <- ns; unlockBinding(nm, ns); assign(nm, f, envir = ns); lockBinding(nm, ns) }
o_curved <- ns$se_row_curved
for (cs in list(list(seed = 2, sg = 1, c0 = 1), list(seed = 2, sg = 1, c0 = 1000), list(seed = 4, sg = 0.1, c0 = 1), list(seed = 4, sg = 0.1, c0 = 3))) {
  set.seed(60000 + cs$seed); n <- 50L
  x <- sort(round(runif(n, 0, 8), 4)); z <- rnorm(n); g <- factor(sample(letters[1:5], n, TRUE))
  y <- cs$c0 * (sin(2 * x) + 0.3 * z + cs$sg * rnorm(n))
  d <- data.frame(x, z, g, y)
  fml <- if (cs$seed == 2) y ~ z + gp(x) else y ~ gp(x, z)
  f <- suppressWarnings(frm(fml, d))
  V <- as.matrix(vcov(f, full = TRUE))
  patch("se_row_curved", function(j, fit, p, u) FALSE)
  f0 <- suppressWarnings(frm(fml, d))
  V0 <- as.matrix(vcov(f0, full = TRUE))
  patch("se_row_curved", o_curved)
  # the closed form's Hessian at the same point
  D2 <- if (cs$seed == 2) outer(x, x, "-")^2 else outer(x, x, "-")^2 + outer(z, z, "-")^2
  X <- if (cs$seed == 2) cbind(1, z) else matrix(1, n, 1)
  nll <- function(p) {
    k <- ncol(X); b <- p[1:k]; K <- exp(2 * p[k + 2]) * exp(-D2 / (2 * exp(2 * p[k + 3]))); diag(K) <- diag(K) * (1 + 1e-12) + exp(2 * p[k + 1])
    R <- chol(K); e <- backsolve(R, y - X %*% b, transpose = TRUE)
    0.5 * sum(e^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
  }
  p0 <- c(f$estimates$beta, f$estimates$betad, f$estimates$theta)
  Hn <- optimHess(p0, nll, control = list(ndeps = rep(1e-3, length(p0))))
  ev <- eigen(Hn, symmetric = TRUE, only.values = TRUE)$values
  sc <- if (all(ev > 0)) sqrt(diag(solve(Hn))) else rep(NA, length(p0))
  fx <- suppressWarnings(fitted(f, newdata = d[1:3, ]))
  cat(sprintf("RSE %s s=%g c=%g conv %d estimates %s\n", deparse(fml), cs$sg, cs$c0, f$opt$convergence, paste(signif(p0, 4), collapse = " ")))
  cat("RSE   shipped SE ", signif(sqrt(diag(V)), 4), "\n")
  cat("RSE   old rule SE", signif(sqrt(diag(V0)), 4), "\n")
  cat("RSE   closed form", signif(sc, 4), " | Hessian eigen", signif(range(ev), 3), "\n")
  flush.console()
}
cat("DONE\n")
