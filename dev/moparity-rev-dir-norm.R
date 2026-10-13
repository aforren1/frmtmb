# Reviewer: the stereographic Dirichlet density integrates to one over
# the chart, and its pushforward to the simplex is the Dirichlet, for
# D = 3 (quadrature) and D = 4 (importance sampling from a Cauchy).
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ld <- function(z, a) {
  frmtmb:::dirichlet_logdens(z, frmtmb:::dirichlet_dist(a, length(a)))
}
# D = 3: polar quadrature on r in (0, inf) via r = tan(t)
a <- c(2, 1.5, 3)
nt <- 600; nphi <- 600
tt <- (seq_len(nt) - 0.5) / nt * (pi / 2)
ph <- (seq_len(nphi) - 0.5) / nphi * 2 * pi
tot <- 0
m1 <- 0
for (t in tt) {
  r <- tan(t)
  jac <- r / cos(t)^2 * (pi / 2 / nt) * (2 * pi / nphi)
  for (p in ph) {
    z <- c(r * cos(p), r * sin(p))
    v <- exp(ld(z, a)) * jac
    tot <- tot + v
    m1 <- m1 + v * frmtmb::mo_simplex(z)[1]
  }
}
cat(sprintf("D = 3, a = (2, 1.5, 3): integral %.6f | E[w1] %.6f (Dirichlet %.6f)\n",
            tot, m1 / tot, a[1] / sum(a)))
# D = 4: importance sampling with a multivariate Cauchy proposal
set.seed(1)
a <- c(0.7, 2, 1, 3)
N <- 2e5
k <- 3
x <- matrix(rnorm(N * k), N)
s <- abs(rnorm(N))
z <- x / s
lq <- lgamma((1 + k) / 2) - lgamma(1 / 2) - (k / 2) * log(pi) -
  ((1 + k) / 2) * log(1 + rowSums(z^2))
lp <- apply(z, 1, ld, a = a)
wt <- exp(lp - lq)
W <- t(apply(z, 1, frmtmb::mo_simplex))
cat(sprintf("D = 4, a = (0.7, 2, 1, 3): integral %.4f (se %.4f) | E[w] %s (Dirichlet %s)\n",
            mean(wt), sd(wt) / sqrt(N),
            paste(format(colSums(W * wt) / sum(wt), digits = 4), collapse = " "),
            paste(format(a / sum(a), digits = 4), collapse = " ")))
