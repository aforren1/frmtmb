# Lane seflags, item 3: the curvature along the smallest direction, at
# the estimate and after Newton steps on the other directions, from
# central differences of the gradient at steps 1e-3, 1e-2 and 1e-1
# (unit-diagonal units), on curved ridges (a + exp(b), its (1 | g)
# variant, c0k) and on identified designs with a small eigenvalue (a
# raw polynomial of degree 5, a near-collinear poisson pair), seeds 1-6.
#   Rscript dev/seflags-ridgecurv.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
quiet <- function(e) suppressMessages(suppressWarnings(e))
probe <- function(f, lab) {
  p <- f$opt$par
  H <- if (!length(f$obj$env$random)) { h <- f$obj$he(p); (h + t(h)) / 2 } else ns$fit_outer_hessian(f)$H
  D <- sqrt(abs(diag(H))); e <- eigen(H / outer(D, D), symmetric = TRUE)
  k <- length(e$values); vk <- e$vectors[, k]; lam <- e$values[k]
  oth <- seq_len(k - 1); U <- e$vectors[, oth, drop = FALSE]; lo <- e$values[oth]
  x <- p; fx <- f$obj$fn(x)
  for (it in 1:30) {
    gx <- as.numeric(f$obj$gr(x)); cc <- drop(crossprod(U, gx / D)); gain <- sum(cc^2/lo)/2
    if (gain < 1e-12) break
    y <- x - drop(U %*% (cc / lo)) / D; fy <- f$obj$fn(y)
    if (!(fy <= fx)) break
    x <- y; fx <- fy
  }
  gt <- sum(as.numeric(f$obj$gr(x)) / D * vk)
  cur <- function(z, eps) { d <- eps * vk / D
    sum((as.numeric(f$obj$gr(z + d)) - as.numeric(f$obj$gr(z - d))) / D * vk) / (2 * eps) }
  cat(sprintf("%-10s lam %.3g | its %d moved %.2g | tangent grad %.2g | curv at p: %s | at p*: %s\n",
    lab, lam, it, max(abs((x - p) * D)), gt,
    paste(signif(sapply(c(1e-3, 1e-2, 1e-1), function(eps) cur(p, eps)), 3), collapse = " "),
    paste(signif(sapply(c(1e-3, 1e-2, 1e-1), function(eps) cur(x, eps)), 3), collapse = " ")))
}
for (s in 1:6) {
  set.seed(s); n <- 60
  dd <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:6, 10)))
  dd$y <- 3 + 0.5 * dd$x + 0.3 * dd$z + rnorm(6, 0, 0.5)[dd$g] + rnorm(n, 0, 0.4)
  probe(quiet(frm(bf(y ~ a + exp(b), a ~ 1 + x, b ~ 1, nl = TRUE), data = dd, start = list(beta = c(1, 0, -1)))), paste("aexpb", s))
  probe(quiet(frm(bf(y ~ a + exp(b), a ~ 1 + x + (1 | g), b ~ 1, nl = TRUE), data = dd, start = list(beta = c(1, 0, -1)))), paste("aexpbre", s))
  set.seed(s); G <- 8
  dn <- data.frame(g = factor(rep(seq_len(G), each = 10))); dn$x <- rnorm(nrow(dn))
  u <- rnorm(G, 0, 0.6)[dn$g]; dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
  f <- tryCatch(quiet(frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE), data = dn)), error = function(e) NULL)
  if (!is.null(f)) probe(f, paste("c0k", s))
  # identified: degree-5 raw polynomial on [1, 2]; poisson near-collinear pair
  set.seed(s); d5 <- data.frame(x = runif(100, 1, 2)); d5$y <- sin(3 * d5$x) + rnorm(100, 0, 0.2)
  probe(quiet(frm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5), data = d5)), paste("poly5", s))
  set.seed(s); x1 <- rnorm(200); dp <- data.frame(x1 = x1, x2 = x1 + rnorm(200, 0, 1e-4)); dp$y <- rpois(200, exp(0.5 + 0.3 * dp$x1))
  probe(quiet(frm(y ~ x1 + x2, family = poisson(), data = dp)), paste("pois", s))
}
