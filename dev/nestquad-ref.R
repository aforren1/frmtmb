# Brute-force references for nested quadrature, written from the model
# alone (plain R densities, no tape, no frmtmb code). Sourced by the
# dev scripts of lane nestquad.

# Gauss-Hermite rule for the standard normal weight (Golub-Welsch),
# written again here so the reference shares no code with the package.
ref_gh <- function(k) {
  i <- seq_len(k - 1)
  J <- matrix(0, k, k)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}

ref_lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }

# log int exp(f(z)) dz by stats::integrate after centering on the mode
# and scaling by the curvature, both found numerically
ref_log_int <- function(f, rel.tol = 1e-11) {
  o <- stats::optimize(function(z) -f(z), c(-60, 60))
  m <- o$minimum
  fm <- f(m)
  h <- 1e-3
  c2 <- -(f(m + h) - 2 * fm + f(m - h)) / h^2
  s <- if (is.finite(c2) && c2 > 0) 1 / sqrt(c2) else 1
  g <- function(x) exp(vapply(m + s * x, f, 1) - fm)
  r <- stats::integrate(g, -40, 40, rel.tol = rel.tol,
                        subdivisions = 10000L, stop.on.error = FALSE)
  if (!r$message %in% "OK") warning("integrate: ", r$message)
  fm + log(s) + log(r$value)
}

# Row log-densities of the families the lane tests, from stats alone.
ref_lpdf <- function(family, y, eta, extra = list()) {
  switch(family,
    bernoulli = stats::dbinom(y, 1, stats::plogis(eta), log = TRUE),
    poisson = stats::dpois(y, exp(eta), log = TRUE),
    gaussian = stats::dnorm(y, eta, extra$sigma, log = TRUE),
    stop("no reference density for ", family))
}

# The marginal log-likelihood of y ~ X beta + (1 | ga) + (1 | ga:gb) by
# nested adaptive integrate(): the integral over each ga effect of the
# product of the integrals over its ga:gb effects. `inner` is the
# ga:gb factor; rows are grouped by `outer` and `inner`.
ref_nested_integrate <- function(family, y, eta0, outer, inner, sd_o,
                                 sd_i, w = rep(1, length(y)),
                                 extra = list(), rel.tol = 1e-10) {
  tot <- 0
  for (a in unique(outer)) {
    ra <- which(outer == a)
    cs <- unique(inner[ra])
    fo <- function(u) {
      s <- stats::dnorm(u, 0, sd_o, log = TRUE)
      for (cc in cs) {
        r <- ra[inner[ra] == cc]
        fi <- function(v) {
          stats::dnorm(v, 0, sd_i, log = TRUE) +
            sum(w[r] * ref_lpdf(family, y[r], eta0[r] + u + v, extra))
        }
        s <- s + ref_log_int(fi, 1e-13)
      }
      s
    }
    tot <- tot + ref_log_int(fo, rel.tol)
  }
  tot
}

# The same integral by non-adaptive Gauss-Hermite on the effects' own
# scale, k nodes per level. Converges as k grows for small clusters.
ref_nested_gh <- function(family, y, eta0, outer, inner, sd_o, sd_i,
                          w = rep(1, length(y)), extra = list(), k = 80) {
  g <- ref_gh(k)
  lw <- log(g$w)
  tot <- 0
  for (a in unique(outer)) {
    ra <- which(outer == a)
    cs <- unique(inner[ra])
    lu <- vapply(sd_o * g$x, function(u) {
      s <- 0
      for (cc in cs) {
        r <- ra[inner[ra] == cc]
        lv <- vapply(sd_i * g$x, function(v) {
          sum(w[r] * ref_lpdf(family, y[r], eta0[r] + u + v, extra))
        }, 1)
        s <- s + ref_lse(lw + lv)
      }
      s
    }, 1)
    tot <- tot + ref_lse(lw + lu)
  }
  tot
}
