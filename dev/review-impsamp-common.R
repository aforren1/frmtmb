# Adversarial review of the impsamp lane: common setup and an independent
# exact marginal, written from the densities and not from the package.
LIB_AFTER <- "C:/Users/adf44/source/r/wt-impsamp-lib"
LIB_BASE <- "C:/Users/adf44/source/r/rellib-r16"
LIB_USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
set_libs <- function(arm = Sys.getenv("REV_ARM", "after")) {
  libs <- switch(arm,
    after = c(LIB_AFTER, LIB_BASE, LIB_USER),
    before = c(LIB_BASE, LIB_USER),
    mutant = c(Sys.getenv("REV_MUTLIB"), LIB_AFTER, LIB_BASE, LIB_USER))
  .libPaths(libs)
  invisible(libs)
}

# Exact log integral over one scalar random effect u ~ N(0, sd^2) of
# prod_i dens(y_i | eta_i + u): an inner mode search then integrate() over
# a window around the mode, so a narrow integrand is not missed.
exact_group_1d <- function(loglik_u, sd) {
  neg <- function(u) -(loglik_u(u) + dnorm(u, 0, sd, log = TRUE))
  o <- optimize(neg, c(-30 * sd - 30, 30 * sd + 30), tol = 1e-8)
  # refine: optimize can miss on a very wide bracket
  o2 <- optim(o$minimum, neg, method = "BFGS",
              control = list(reltol = 1e-12))
  mode <- if (o2$value < o$objective) o2$par else o$minimum
  v0 <- min(o$objective, o2$value)
  h <- 1e-3
  cur <- (neg(mode + h) - 2 * neg(mode) + neg(mode - h)) / h^2
  s <- if (is.finite(cur) && cur > 0) 1 / sqrt(cur) else sd
  f <- function(x) vapply(x, function(xi) exp(v0 - neg(mode + s * xi)), 0)
  val <- integrate(f, -Inf, Inf, rel.tol = 1e-10, abs.tol = 0,
                   subdivisions = 2000L)$value
  log(val) + log(s) - v0
}

# Exact log integral over a 2-d block u ~ N(0, Sigma): product Gauss-Hermite
# of order m around the mode with the mode's Hessian. Not the package's.
gh_rule <- function(m) {
  j <- seq_len(m - 1L)
  b <- sqrt(j / 2)
  J <- matrix(0, m, m)
  J[cbind(j, j + 1L)] <- b
  J[cbind(j + 1L, j)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1L, ]^2)
}
exact_group_2d <- function(loglik_u, Sigma, m = 40L) {
  Pinv <- solve(Sigma)
  ld <- as.numeric(determinant(Sigma)$modulus)
  neg <- function(u) {
    -(loglik_u(u) - 0.5 * sum(u * (Pinv %*% u)) - log(2 * pi) - 0.5 * ld)
  }
  o <- optim(c(0, 0), neg, method = "BFGS", hessian = FALSE,
             control = list(reltol = 1e-12, maxit = 500))
  o <- optim(o$par, neg, method = "BFGS", hessian = TRUE,
             control = list(reltol = 1e-14, maxit = 500))
  L <- t(chol(solve(o$hessian)))
  g <- gh_rule(m)
  ng <- expand.grid(i = seq_len(m), j = seq_len(m))
  vals <- mapply(function(i, j) {
    z <- c(g$x[i], g$x[j])
    u <- o$par + as.vector(sqrt(2) * L %*% z)
    log(g$w[i]) + log(g$w[j]) + sum(z^2) - neg(u)
  }, ng$i, ng$j)
  mx <- max(vals)
  mx + log(sum(exp(vals - mx))) + log(2) + sum(log(diag(L)))
}

log_mean_exp <- function(a) {
  m <- max(a)
  m + log(mean(exp(a - m)))
}
