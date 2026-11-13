# Review of lane arspat: the Laplace error of fcor()/sar() + mi() on
# student() at a moderate nu. The objective (Laplace over the missing
# values) is read at the fit's estimate with nu moved to 3, 5 and 10,
# against the exact multivariate-t margin over the observed rows.
source("dev/reviews/arspat-review-common.R")
suppressPackageStartupMessages(library(frmtmb.spatial))

sar_W <- function(n, seed = 5, k = 4) {
  set.seed(seed)
  xy <- cbind(runif(n), runif(n))
  D <- as.matrix(dist(xy)); diag(D) <- Inf
  K <- t(apply(D, 1, function(r) as.numeric(rank(r) <= k)))
  K <- pmax(K, t(K))
  K / rowSums(K)
}

run <- function(kind, nmiss, n = 60, seed = 3) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  if (kind == "fcor") {
    M <- 0.5^abs(outer(1:n, 1:n, "-"))
    Sbase <- M
    d$y <- 1 + 0.5 * d$x + as.vector(t(chol(M)) %*% rnorm(n)) * 0.8 * 1.6
    fo <- y | mi() ~ x + fcor(M)
    d2 <- list(M = M)
  } else {
    W <- sar_W(n)
    A0 <- diag(n) - 0.4 * W
    d$y <- 1 + 0.5 * d$x + as.vector(solve(A0, rnorm(n, 0, 0.7))) * 1.6
    fo <- y | mi() ~ x + sar(W, type = "error")
    d2 <- list(W = W)
  }
  miss <- sort(sample(n, nmiss))
  d$y[miss] <- NA
  f <- frm(fo, data = d, data2 = d2, family = student())
  lab <- ns$outer_par_labels(f)
  par <- f$opt$par
  inu <- grep("^nu", lab)
  ok <- !is.na(d$y)
  for (nu in c(3, 5, 10, 30)) {
    p <- par
    p[inu] <- log(nu - 1)
    lap <- -f$obj$fn(p)
    b <- p[grep("^b|Intercept|^x$", lab)][1:2]
    beta <- unname(p[which(ns$outer_par_map(f)$comp == "beta")])
    sig <- exp(unname(p[grep("^sigma", lab)]))
    mu <- beta[1] + beta[2] * d$x
    if (kind == "fcor") {
      S <- sig^2 * Sbase
    } else {
      ac <- f$frame$autocor[[1L]]
      th <- unname(p[which(ns$outer_par_map(f)$comp == "thetaac")])
      rho <- unname(ns$autocor_natural(th, ac))
      A <- diag(n) - rho * W
      S <- sig^2 * solve(crossprod(A))
    }
    ex <- mvtnorm::dmvt(d$y[ok], mu[ok], S[ok, ok], df = nu, log = TRUE)
    cat(sprintf("%-4s nmiss %2d nu %4.0f: laplace %.6f exact %.6f diff %+.3e\n",
                kind, nmiss, nu, lap, ex, lap - ex))
  }
  invisible(lab)
}
chk("labels", print(run("fcor", 3)))
for (k in c(10, 20)) chk(paste("fcor", k), run("fcor", k))
for (k in c(3, 10, 20)) chk(paste("sar", k), run("sar", k))
