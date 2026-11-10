# Part C: can the flat direction of the ridge be separated from the weak
# c0 direction beside it? Rayleigh-Ritz on the two smallest directions of
# the retaken Hessian, with the Hessian times each direction measured by
# central differences of the gradient after the larger directions are
# polished (the procedure of se_ridge_curvature()), at steps h and 10 h
# extrapolated to zero.
#   CODE8_ARM=before Rscript dev/code8-ridge-rr.R <h> <j...>
args <- commandArgs(TRUE)
h0 <- as.numeric(args[1]); js <- as.integer(args[-1])
src <- readLines("dev/code8-ridge-exp.R")
a <- c("retake", "rr")
eval(parse(text = src[grep("^source", src):grep("^set.seed", src)[1] - 1]))
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
y0 <- dn$yn
form <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
for (j in js) {
  dn$yn <- y0 * (1 + 2e-15 * j)
  f <- quiet(frm(form, data = dn))
  p <- f$opt$par
  H <- fns$fit_outer_hessian(f)$H
  D <- sqrt(abs(diag(H)))
  e <- eigen(H / outer(D, D), symmetric = TRUE)
  np <- length(p)
  S <- c(np - 1, np)
  U <- e$vectors[, -S]; lo <- e$values[-S]
  g <- function(x) drop(f$obj$gr(x)) / D
  x <- p
  gx <- g(x)
  for (it in 1:20) {
    y <- x - drop(U %*% (drop(crossprod(U, gx)) / lo)) / D
    gy <- g(y)
    if (!(sum(crossprod(U, gy)^2 / lo) < sum(crossprod(U, gx)^2 / lo))) break
    x <- y; gx <- gy
  }
  Hv <- function(v, h) {
    (g(x + h * v / D) - g(x - h * v / D)) / (2 * h)
  }
  V <- e$vectors[, S]
  W1 <- sapply(1:2, function(i) Hv(V[, i], h0))
  W2 <- sapply(1:2, function(i) Hv(V[, i], 10 * h0))
  W <- (100 * W1 - W2) / 99
  B <- crossprod(V, W)
  cat(sprintf("== j %d  polished moved %.3g\n", j, sqrt(sum(((x - p) * D)^2))))
  cat(" B:", format(B, digits = 4), " asym", format(abs(B[1, 2] - B[2, 1]),
                                                     digits = 3), "\n")
  eb <- eigen((B + t(B)) / 2, symmetric = TRUE)
  Wv <- V %*% eb$vectors
  v <- c(0, 0, p[3], -p[4], 0, 1); tu <- D * v; tu <- tu / sqrt(sum(tu^2))
  cat(" eig:", format(eb$values, digits = 4), "\n flat vec:",
      format(Wv[, 2], digits = 4), "\n true tangent:", format(tu, digits = 4),
      "\n Hessian's own:", format(e$vectors[, np], digits = 4), "\n")
}
cat("DONE\n")
