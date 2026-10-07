# Lane surface, punch round 1: the accuracy of the n-point Gauss-Hermite
# SD of plogis(Z), Z ~ N(mu, s^2), against integrate(), over a grid.
#   Rscript dev/surface-p1-gh.R > dev/surface-out/p1-gh.txt
source("dev/surface-env.R"); surface_env("lane")
gh <- function(n) {
  J <- matrix(0, n, n); off <- sqrt(seq_len(n - 1L))
  J[cbind(1:(n - 1), 2:n)] <- off; J[cbind(2:n, 1:(n - 1))] <- off
  e <- eigen(J, symmetric = TRUE); list(x = e$values, w = e$vectors[1, ]^2)
}
sd_gh <- function(mu, s, g) {
  v <- plogis(mu + s * g$x); m <- sum(g$w * v); sqrt(sum(g$w * (v - m)^2))
}
sd_ref <- function(mu, s) {
  f <- function(z) plogis(z) * dnorm(z, mu, s)
  m <- integrate(f, mu - 14 * s, mu + 14 * s, rel.tol = 1e-13,
                 subdivisions = 1000L)$value
  sqrt(integrate(function(z) (plogis(z) - m)^2 * dnorm(z, mu, s),
                 mu - 14 * s, mu + 14 * s, rel.tol = 1e-13,
                 subdivisions = 1000L)$value)
}
grid <- expand.grid(mu = c(-20, -16, -10, -5, -2, 0), s = c(0.1, 0.5, 1, 2, 4, 6, 10))
for (n in c(40, 100, 200)) {
  g <- gh(n)
  r <- mapply(function(mu, s) sd_gh(mu, s, g) / sd_ref(mu, s), grid$mu, grid$s)
  cat(sprintf("n = %3d: max |ratio - 1| %.3g; at s <= 4: %.3g\n", n,
              max(abs(r - 1)), max(abs(r[grid$s <= 4] - 1))))
}
