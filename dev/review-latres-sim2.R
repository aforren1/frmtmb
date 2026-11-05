# Attack 2 and 3: dispersion families (Gamma, negbinomial) and the cov = FALSE
# forms, against the exact marginal by forward recursion on a grid.
#
# cov = TRUE:  err is AR(1); state err_t, one-dimensional as in review-latres-sim.R.
# cov = FALSE: mu_t = g(eta_t + s z_t) + ar s z_{t-1}  (on_response)
#              eta_t + s z_t + ar s z_{t-1}            (link scale)
#              z iid; state z_t. p(y_t | z_t, z_{t-1}) is a matrix on the grid.
#              A mean outside the support has density 0, as Stan rejects it.
# Usage: Rscript review-latres-sim2.R <design> <seed_from> <seed_to> <outfile>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
args <- commandArgs(TRUE)
design <- args[[1]]
s1 <- as.integer(args[[2]])
s2 <- as.integer(args[[3]])
outf <- args[[4]]
designs <- list(
  gamma_T = list(kind = "T", fam = Gamma(link = "log"), dist = "gamma", shape = 4, b0 = 0.5, b1 = 0.4, s = 0.5, ar = 0.5),
  nb_T = list(kind = "T", fam = negbinomial(), dist = "nb", shape = 4, b0 = 1.5, b1 = 0.4, s = 0.5, ar = 0.5),
  gamma_F = list(kind = "F", fam = Gamma(link = "log"), dist = "gamma", post = TRUE, inv = exp, shape = 4, b0 = 0.5, b1 = 0.4, s = 0.5, ar = 0.5),
  probit_F = list(kind = "F", fam = bernoulli("probit"), dist = "bern", post = TRUE, inv = pnorm, b0 = 0.2, b1 = 0.5, s = 0.5, ar = 0.2),
  pois_F = list(kind = "F", fam = poisson(), dist = "pois", post = FALSE, inv = exp, b0 = 0.5, b1 = 0.4, s = 0.5, ar = 0.5),
  bern_F = list(kind = "F", fam = bernoulli(), dist = "bern", post = FALSE, inv = plogis, b0 = 0.2, b1 = 0.5, s = 1, ar = 0.5))
D <- designs[[design]]
G <- 25
K <- 6
ngrid <- if (D$kind == "F") 61 else 121

dens <- function(y, mu, shape) {
  v <- switch(D$dist,
    gamma = dgamma(y, shape, shape / mu),
    nb = dnbinom(y, size = shape, mu = mu),
    pois = dpois(y, mu),
    bern = dbinom(y, 1, mu))
  v[!is.finite(v)] <- 0
  v
}
simd <- function(seed) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = seq_len(G))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- rnorm(n)
  d$y <- NA_real_
  for (g in seq_len(G)) {
    idx <- which(d$g == g)
    if (D$kind == "T") {
      e <- as.numeric(arima.sim(list(ar = D$ar), K, sd = D$s))
      mu <- D$fam$linkinv(D$b0 + D$b1 * d$x[idx] + e)
    } else {
      repeat {
        z <- rnorm(K)
        lag <- c(0, z[-K])
        mu <- if (D$post) {
          D$inv(D$b0 + D$b1 * d$x[idx] + D$s * z) + D$ar * D$s * lag
        } else {
          D$inv(D$b0 + D$b1 * d$x[idx] + D$s * z + D$ar * D$s * lag)
        }
        ok <- switch(D$dist, bern = all(mu > 0 & mu < 1), all(mu > 0))
        if (ok) break
      }
    }
    d$y[idx] <- switch(D$dist,
      gamma = rgamma(K, D$shape, D$shape / mu),
      nb = rnbinom(K, size = D$shape, mu = mu),
      pois = rpois(K, mu),
      bern = rbinom(K, 1, mu))
  }
  d$g <- factor(d$g)
  d
}
inv_link <- D$fam$linkinv
exact_T <- function(par, d, ng = ngrid) {
  b0 <- par[1]; b1 <- par[2]; s <- exp(par[3]); th <- par[4]
  shape <- if (length(par) > 4) exp(par[5]) else NA
  ar <- th / sqrt(1 + th^2)
  msd <- s * sqrt(1 + th^2)
  u <- seq(-6.5, 6.5, length.out = ng)
  h <- u[2] - u[1]
  grid <- u * msd
  Tm <- outer(grid, grid, function(a, b) dnorm(b, ar * a, s)) * (h * msd)
  ll <- 0
  for (idx in split(seq_len(nrow(d)), d$g)) {
    f <- dnorm(u) * h
    for (k in seq_along(idx)) {
      i <- idx[k]
      if (k > 1L) f <- as.vector(crossprod(f, Tm))
      f <- f * dens(d$y[i], inv_link(b0 + b1 * d$x[i] + grid), shape)
      c0 <- sum(f)
      if (!is.finite(c0) || c0 <= 0) return(1e10)
      ll <- ll + log(c0)
      f <- f / c0
    }
  }
  -ll
}
exact_F <- function(par, d, ng = ngrid) {
  b0 <- par[1]; b1 <- par[2]; s <- exp(par[3]); ar <- 2 / (1 + exp(-par[4])) - 1
  shape <- if (length(par) > 4) exp(par[5]) else NA
  u <- seq(-6, 6, length.out = ng)
  h <- u[2] - u[1]
  pz <- dnorm(u) * h
  ll <- 0
  lagm <- ar * s * u
  for (idx in split(seq_len(nrow(d)), d$g)) {
    eta <- b0 + b1 * d$x[idx[1]] + s * u
    mu <- if (D$post) D$inv(eta) else D$inv(eta)
    f <- pz * dens(d$y[idx[1]], mu, shape)
    c0 <- sum(f)
    if (!is.finite(c0) || c0 <= 0) return(1e10)
    ll <- ll + log(c0)
    f <- f / c0
    for (i in idx[-1]) {
      eta <- b0 + b1 * d$x[i] + s * u       # by z_t
      M <- if (D$post) {
        outer(lagm, D$inv(eta), "+")        # z_{t-1} in rows, z_t in columns
      } else {
        D$inv(outer(lagm, eta, "+"))
      }
      P <- matrix(dens(d$y[i], as.vector(M), shape), ng, ng)
      f <- as.vector(crossprod(f, P)) * pz
      c0 <- sum(f)
      if (!is.finite(c0) || c0 <= 0) return(1e10)
      ll <- ll + log(c0)
      f <- f / c0
    }
  }
  -ll
}
exact_nll <- if (D$kind == "T") exact_T else exact_F

one <- function(seed) {
  d <- simd(seed)
  form <- if (D$kind == "T") y ~ x + ar(t, g, cov = TRUE) else y ~ x + ar(t, g)
  t0 <- proc.time()[[3]]
  m <- msg_of(frm(form, data = d, family = D$fam))
  tfit <- proc.time()[[3]] - t0
  has_shape <- !is.null(D$shape)
  out <- list(design = design, seed = seed, tfit = tfit, err = substr(m$err %||% "", 1, 100))
  st <- c(D$b0, D$b1, log(D$s),
          if (D$kind == "T") D$ar / sqrt(1 - D$ar^2) else qlogis((D$ar + 1) / 2),
          if (has_shape) log(D$shape))
  frm_par <- rep(NA_real_, length(st))
  if (is.null(m$err)) {
    f <- m$val
    th <- f$estimates$theta
    frm_par <- c(f$estimates$beta[1:2], th[1], th[2],
                 if (has_shape) f$opt$par[names(f$opt$par) == "betad"])
    out$ll <- as.numeric(logLik(f))
    out$code <- f$opt$convergence
    vc <- tryCatch(as.matrix(vcov(f, full = TRUE)), error = function(e) NULL)
    out$se_b1 <- tryCatch(sqrt(vc[2, 2]), error = function(e) NA)
    out$se_ls <- tryCatch(sqrt(vc["theta_1", "theta_1"]), error = function(e) NA)
    out$nwarn <- length(m$warn)
    out$warn <- substr(paste(m$warn, collapse = " | "), 1, 120)
    out$frm_b0 <- frm_par[1]; out$frm_b1 <- frm_par[2]; out$frm_ls <- frm_par[3]; out$frm_t <- frm_par[4]
    out$frm_ar <- if (D$kind == "T") frm_par[4] / sqrt(1 + frm_par[4]^2) else 2 / (1 + exp(-frm_par[4])) - 1
    if (has_shape) out$frm_lshape <- frm_par[5]
    pf <- frm_par
    if (D$kind == "T") pf[4] <- pmin(pmax(pf[4], -8), 8)
    if (has_shape) pf[5] <- pmin(pf[5], log(1e6))
    out$x_at_frm <- -exact_nll(pf, d)
  }
  lo <- c(-Inf, -Inf, log(0.01), if (D$kind == "T") -8 else -12, if (has_shape) log(0.05))
  hi <- c(Inf, Inf, log(10), if (D$kind == "T") 8 else 12, if (has_shape) log(1e4))
  t0 <- proc.time()[[3]]
  ex <- tryCatch(nlminb(st, exact_nll, d = d, lower = lo, upper = hi,
                        control = list(iter.max = 100, eval.max = 300)),
                 error = function(e) NULL)
  if (!is.null(ex)) {
    out$tex <- proc.time()[[3]] - t0
    out$x_b0 <- ex$par[1]; out$x_b1 <- ex$par[2]; out$x_ls <- ex$par[3]; out$x_t <- ex$par[4]
    out$x_ar <- if (D$kind == "T") ex$par[4] / sqrt(1 + ex$par[4]^2) else 2 / (1 + exp(-ex$par[4])) - 1
    if (has_shape) out$x_lshape <- ex$par[5]
    out$x_ll <- -ex$objective
    out$x_code <- ex$convergence
  }
  out
}
cols <- c("design","seed","tfit","err","ll","code","se_b1","se_ls","nwarn","warn","frm_b0","frm_b1","frm_ls","frm_t","frm_ar","frm_lshape","x_at_frm","tex","x_b0","x_b1","x_ls","x_t","x_ar","x_lshape","x_ll","x_code")
first <- TRUE
for (sd in s1:s2) {
  o <- one(sd)
  for (cn in setdiff(cols, names(o))) o[[cn]] <- NA
  write.table(as.data.frame(o[cols]), outf, sep = "\t", row.names = FALSE, quote = TRUE, col.names = first, append = !first)
  first <- FALSE
}
