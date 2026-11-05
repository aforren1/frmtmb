# The grid's cells, its data, and an adaptive Gauss-Hermite reference
# that shares no code with frmtmb: R's own density functions, a per-group
# Newton solve in R, and tensor Gauss-Hermite nodes. It is here to judge
# frmtmb, so nothing in it reads the package's internals.

impsamp_cells <- function() {
  expand.grid(
    family = c("binomial", "poisson", "bernoulli", "negbinomial", "gaussian",
               "Gamma"),
    shape = c("int", "slope", "twoblock"),
    size = c("10x20", "40x5", "30x3"),
    seed = 1:3,
    N = c(64L, 200L, 512L, 2000L),
    stringsAsFactors = FALSE)
}

impsamp_size <- function(size) {
  as.integer(strsplit(size, "x", fixed = TRUE)[[1L]])
}

impsamp_data <- function(family, shape, size, seed) {
  sz <- impsamp_size(size)
  G <- sz[[1L]]
  m <- sz[[2L]]
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (!is.null(old)) assign(".Random.seed", old, globalenv()))
  # the cell's own stream, so a data set does not depend on which chunk
  # or arm drew it
  set.seed(1000L * seed + match(family, impsamp_cells()$family[1:6]) +
             100L * match(size, c("10x20", "40x5", "30x3")))
  g <- factor(rep(seq_len(G), each = m))
  x <- round(rnorm(G * m), 3)
  sd1 <- 0.8
  sd2 <- 0.4
  rho <- if (shape == "slope") 0.3 else 0
  a <- rnorm(G)
  b <- rho * a + sqrt(1 - rho^2) * rnorm(G)
  u1 <- sd1 * a
  u2 <- if (shape == "int") rep(0, G) else sd2 * b
  eta <- 0.3 * x + u1[g] + u2[g] * x
  n <- G * m
  y <- switch(family,
    binomial = rbinom(n, 10L, plogis(eta)),
    bernoulli = rbinom(n, 1L, plogis(eta - 0.5)),
    poisson = rpois(n, exp(eta + 0.5)),
    negbinomial = rnbinom(n, mu = exp(eta + 1), size = 2),
    gaussian = eta + 0.2 + rnorm(n, 0, 0.5),
    Gamma = rgamma(n, shape = 3, rate = 3 / exp(eta + 0.5)))
  data.frame(y = y, x = x, g = g, ntr = 10L)
}

impsamp_formula <- function(family, shape) {
  re <- switch(shape, int = "(1 | g)", slope = "(1 + x | g)",
               twoblock = "(1 | g) + (0 + x | g)")
  lhs <- if (family == "binomial") "y | trials(ntr)" else "y"
  frmtmb::bf(stats::as.formula(paste(lhs, "~ x +", re)))
}

impsamp_family <- function(family) {
  switch(family, binomial = stats::binomial(), bernoulli = frmtmb::bernoulli(),
         poisson = stats::poisson(), negbinomial = frmtmb::negbinomial(),
         gaussian = stats::gaussian(), Gamma = stats::Gamma(link = "log"))
}

# the row log density at eta, with one log dispersion
agh_logdens <- function(family, y, eta, ntr, ld) {
  switch(family,
    binomial = dbinom(y, ntr, plogis(eta), log = TRUE),
    bernoulli = dbinom(y, 1L, plogis(eta), log = TRUE),
    poisson = dpois(y, exp(eta), log = TRUE),
    negbinomial = dnbinom(y, mu = exp(eta), size = exp(ld), log = TRUE),
    gaussian = dnorm(y, eta, exp(ld), log = TRUE),
    Gamma = dgamma(y, shape = exp(ld), rate = exp(ld) / exp(eta), log = TRUE))
}
agh_has_disp <- function(family) family %in% c("negbinomial", "gaussian", "Gamma")

# group-level random design: Zg %*% u
agh_Sigma <- function(shape, th) {
  if (shape == "int") return(matrix(exp(2 * th[[1L]]), 1, 1))
  s1 <- exp(th[[1L]])
  s2 <- exp(th[[2L]])
  r <- if (shape == "slope") tanh(th[[3L]]) else 0
  matrix(c(s1^2, r * s1 * s2, r * s1 * s2, s2^2), 2, 2)
}

agh_nodes <- function(k) {
  # Golub-Welsch for the physicists' weight exp(-x^2)
  a <- sqrt(seq_len(k - 1L) / 2)
  J <- matrix(0, k, k)
  J[cbind(seq_len(k - 1L), 2:k)] <- a
  J[cbind(2:k, seq_len(k - 1L))] <- a
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1L, ]^2)
}

# one group's log marginal likelihood
agh_group <- function(family, y, Zg, etaf, ntr, ld, Sigma, nodes) {
  q <- ncol(Zg)
  Pinv <- solve(Sigma)
  f <- function(u) {
    -sum(agh_logdens(family, y, etaf + as.vector(Zg %*% u), ntr, ld)) +
      0.5 * sum(u * (Pinv %*% u)) +
      0.5 * as.numeric(determinant(Sigma, TRUE)$modulus) +
      0.5 * q * log(2 * pi)
  }
  # the mode, by damped Newton on a numeric gradient of an analytic-free
  # f: this is a reference, so it stays simple and slow rather than clever
  u <- optim(rep(0, q), f, method = "BFGS",
             control = list(reltol = 1e-12, maxit = 500))$par
  h <- 1e-4
  H <- matrix(0, q, q)
  for (i in seq_len(q)) for (j in seq_len(i)) {
    ei <- replace(numeric(q), i, h)
    ej <- replace(numeric(q), j, h)
    H[i, j] <- H[j, i] <- (f(u + ei + ej) - f(u + ei - ej) - f(u - ei + ej) +
                             f(u - ei - ej)) / (4 * h^2)
  }
  R <- tryCatch(chol(H), error = function(e) chol(H + diag(1e-6, q)))
  L <- backsolve(R, diag(q))
  # the tensor grid, in log space
  if (q == 1L) {
    z <- matrix(nodes$x, ncol = 1L)
    lw <- log(nodes$w)
  } else {
    gr <- expand.grid(i = seq_along(nodes$x), j = seq_along(nodes$x))
    z <- cbind(nodes$x[gr$i], nodes$x[gr$j])
    lw <- log(nodes$w[gr$i]) + log(nodes$w[gr$j])
  }
  vals <- vapply(seq_len(nrow(z)), function(i) {
    ui <- u + sqrt(2) * as.vector(L %*% z[i, ])
    -f(ui) + sum(z[i, ]^2) + lw[[i]]
  }, 0)
  mx <- max(vals)
  mx + log(sum(exp(vals - mx))) + 0.5 * q * log(2) + sum(log(diag(L)))
}

agh_ll <- function(family, shape, d, par, k) {
  # par = beta0, beta1, theta..., [log dispersion]
  nth <- switch(shape, int = 1L, slope = 3L, twoblock = 2L)
  beta <- par[1:2]
  th <- par[2L + seq_len(nth)]
  ld <- if (agh_has_disp(family)) par[[3L + nth]] else 0
  Sigma <- agh_Sigma(shape, th)
  nodes <- agh_nodes(k)
  etaf <- beta[[1L]] + beta[[2L]] * d$x
  ll <- 0
  for (lv in levels(d$g)) {
    w <- d$g == lv
    Zg <- if (shape == "int") matrix(1, sum(w), 1) else cbind(1, d$x[w])
    ll <- ll + agh_group(family, d$y[w], Zg, etaf[w], d$ntr[w], ld, Sigma,
                         nodes)
  }
  ll
}

impsamp_start <- function(family, shape, lap) {
  cf <- lap$opt$par
  b <- unname(cf[names(cf) == "beta"])
  # (see impsamp_sds)
  # from the printed scale rather than the chart, so the reference never
  # reads the package's internal coordinates; a singular start is pulled
  # back because the reference's own chart cannot sit on the boundary
  sds <- impsamp_sds(lap)
  th <- switch(shape, int = log(max(sds[[1L]], 0.02)),
    slope = c(log(max(sds[[1L]], 0.02)), log(max(sds[[2L]], 0.02)),
              atanh(max(min(sds[[3L]], 0.99999), -0.99999))),
    twoblock = log(pmax(sds[1:2], 0.02)))
  c(b, th, if (agh_has_disp(family)) impsamp_disp_start(family, lap))
}

# sd, sd, cor from VarCorr()'s documented table
impsamp_sds <- function(fit) {
  vc <- frmtmb::VarCorr(fit)[["g"]]
  s <- unname(vc[["sd"]][, "Estimate"])
  co <- vc[["cor"]]
  c(s, if (!is.null(co)) co[2L, "Estimate", 1L] else 0)
}

impsamp_disp_start <- function(family, fit) {
  # betad is the dispersion's one coefficient under its log link (shape
  # for Gamma and negbinomial, sigma for gaussian)
  log_disp <- unname(fit$opt$par[names(fit$opt$par) == "betad"])
  log_disp[[1L]]
}

# The exact maximum, by optim() on the 2-node-per-unit reference
impsamp_ref <- function(family, shape, d, lap) {
  k <- if (shape == "int") 25L else 11L
  st <- impsamp_start(family, shape, lap)
  fn <- function(p) {
    v <- tryCatch(-agh_ll(family, shape, d, p, k), error = function(e) 1e10)
    if (is.finite(v)) v else 1e10
  }
  o <- optim(st, fn, method = "BFGS",
             control = list(reltol = 1e-9, maxit = 200))
  list(ll = -o$value, par = o$par)
}

# The exact marginal likelihood at a fit's own estimates, so a reported
# logLik is judged against the truth at its own point
impsamp_exact_at <- function(family, shape, d, fit) {
  k <- if (shape == "int") 25L else 11L
  agh_ll(family, shape, d, impsamp_start(family, shape, fit), k)
}
