# M4: which post-link cov = FALSE forms, and which guarded link-scale forms, does
# Laplace read true? Against the exact marginal by forward recursion on a grid
# (dev/review-latres-sim2.R's construction, its 61-point chain for cov = FALSE and
# 121 for cov = TRUE). Also the objective's dependence on the inner start, after
# dev/review-latres-hyst.R, at the fit's outer point.
#
# Usage: Rscript dev/latres-m4.R <lib> <design> <seed_from> <seed_to> <outfile>
# Data seed = seed. G = 25 groups of K = 6.
args <- commandArgs(TRUE)
.libPaths(c(args[[1]], "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
design <- args[[2]]; s1 <- as.integer(args[[3]]); s2 <- as.integer(args[[4]])
outf <- args[[5]]
`%||%` <- function(a, b) if (is.null(a)) b else a
gam <- function(l) Gamma(link = l)
dd <- list(
  # the post-link forms, kind F
  gamma_log = list("F", gam("log"), "gamma", TRUE, 0.5, 0.4, 0.5, 0.5, shape = 4),
  gamma_inv = list("F", gam("inverse"), "gamma", TRUE, 0.5, 0.05, 0.05, 0.3, shape = 4),
  beta_logit = list("F", Beta(), "beta", TRUE, 0.2, 0.3, 0.4, 0.3, shape = 10),
  weibull_log = list("F", weibull(), "weibull", TRUE, 0.8, 0.3, 0.4, 0.4, shape = 2),
  expo_log = list("F", exponential(), "exp", TRUE, 0.8, 0.3, 0.4, 0.4),
  invg_log = list("F", inverse.gaussian(link = "log"), "invg", TRUE, 0.5, 0.3, 0.4, 0.4, shape = 4),
  probit = list("F", bernoulli("probit"), "bern", TRUE, 0.2, 0.5, 0.5, 0.2),
  cloglog = list("F", bernoulli("cloglog"), "bern", TRUE, -0.5, 0.5, 0.5, 0.2),
  binom_probit = list("F", binomial("probit"), "binom", TRUE, 0.2, 0.5, 0.5, 0.2, tr = 10),
  pois_sqrt = list("F", poisson("sqrt"), "pois", TRUE, 2, 0.3, 0.3, 0.3),
  pois_sqrt_lo = list("F", poisson("sqrt"), "pois", TRUE, 1, 0.2, 0.4, 0.4),
  pois_sqrt_big = list("F", poisson("sqrt"), "pois", TRUE, 3, 0.3, 0.6, 0.5),
  pois_cens = list("F", poisson(), "pois", TRUE, 1, 0.4, 0.5, 0.5, cens = TRUE),
  gamma_big = list("F", gam("log"), "gamma", TRUE, 3, 0.3, 0.3, 0.3, shape = 20),
  # link-scale controls, kind F
  pois_ctl = list("F", poisson(), "pois", FALSE, 0.5, 0.4, 0.5, 0.5),
  # the guarded link-scale forms, kind T: a latent can cross the support
  pois_id = list("T", poisson("identity"), "pois", NA, 8, 1, 1, 0.5),
  gamma_id = list("T", gam("identity"), "gamma", NA, 8, 1, 1, 0.5, shape = 4),
  gamma_invT = list("T", gam("inverse"), "gamma", NA, 0.5, 0.05, 0.05, 0.5, shape = 4),
  binom_id = list("T", binomial("identity"), "binom", NA, 0.5, 0.05, 0.05, 0.5, tr = 10),
  pois_idlo = list("T", poisson("identity"), "pois", NA, 2, 0.5, 1, 0.5))
Z <- dd[[design]]
kind <- Z[[1]]; fam <- Z[[2]]; dist <- Z[[3]]; post <- Z[[4]]
b0 <- Z[[5]]; b1 <- Z[[6]]; s0 <- Z[[7]]; ar0 <- Z[[8]]
shape0 <- Z[["shape"]]; tr <- Z[["tr"]] %||% 1L; cens <- isTRUE(Z[["cens"]])
G <- 25L; K <- 6L
ngrid <- if (kind == "F") 61L else 121L
inv <- fam$linkinv

inside <- function(mu) {
  switch(dist, bern = , binom = , beta = mu > 0 & mu < 1, mu > 0)
}
dens <- function(y, mu, ex, cen = 0L) {
  v <- rep(0, length(mu))
  ok <- is.finite(mu) & inside(mu)
  m <- mu[ok]
  sh <- if (is.null(ex) || is.na(ex)) NA else exp(ex)
  v[ok] <- if (cen == 1L) {
    ppois(y, m, lower.tail = FALSE)
  } else switch(dist,
    gamma = dgamma(y, sh, sh / m),
    beta = dbeta(y, m * sh, (1 - m) * sh),
    weibull = dweibull(y, sh, m / gamma(1 + 1 / sh)),
    exp = dexp(y, 1 / m),
    invg = sqrt(sh / (2 * pi * y^3)) * exp(-sh * (y - m)^2 / (2 * m^2 * y)),
    pois = dpois(y, m),
    bern = dbinom(y, 1, m),
    binom = dbinom(y, tr, m))
  v[!is.finite(v)] <- 0
  v
}
rdist <- function(n, mu) {
  switch(dist,
    gamma = rgamma(n, shape0, shape0 / mu),
    beta = pmin(pmax(rbeta(n, mu * shape0, (1 - mu) * shape0), 1e-6), 1 - 1e-6),
    weibull = rweibull(n, shape0, mu / gamma(1 + 1 / shape0)),
    exp = rexp(n, 1 / mu),
    invg = {
      nu <- rnorm(n)^2
      x <- mu + mu^2 * nu / (2 * shape0) -
        mu / (2 * shape0) * sqrt(4 * mu * shape0 * nu + mu^2 * nu^2)
      ifelse(runif(n) <= mu / (mu + x), x, mu^2 / x)
    },
    pois = rpois(n, mu),
    bern = rbinom(n, 1, mu),
    binom = rbinom(n, tr, mu))
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
    repeat {
      if (kind == "T") {
        e <- as.numeric(arima.sim(list(ar = ar0), K, sd = s0))
        mu <- inv(b0 + b1 * d$x[idx] + e)
      } else {
        z <- rnorm(K)
        lag <- c(0, z[-K])
        mu <- if (post) {
          inv(b0 + b1 * d$x[idx] + s0 * z) + ar0 * s0 * lag
        } else {
          inv(b0 + b1 * d$x[idx] + s0 * z + ar0 * s0 * lag)
        }
      }
      if (all(inside(mu))) break
    }
    d$y[idx] <- rdist(K, mu)
  }
  d$cen <- 0L
  if (cens) {
    cut <- 6
    d$cen <- as.integer(d$y > cut)
    d$y[d$cen == 1L] <- cut
  }
  if (dist == "binom") d$tr <- tr
  d$g <- factor(d$g)
  d
}
has_ex <- !is.null(shape0)
exact_T <- function(par, d, ng = ngrid) {
  bb0 <- par[1]; bb1 <- par[2]; s <- exp(par[3]); th <- par[4]
  ex <- if (has_ex) par[5] else NA
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
      f <- f * dens(d$y[i], inv(bb0 + bb1 * d$x[i] + grid), ex)
      c0 <- sum(f)
      if (!is.finite(c0) || c0 <= 0) return(1e10)
      ll <- ll + log(c0)
      f <- f / c0
    }
  }
  -ll
}
exact_F <- function(par, d, ng = ngrid) {
  bb0 <- par[1]; bb1 <- par[2]; s <- exp(par[3]); ar <- 2 / (1 + exp(-par[4])) - 1
  ex <- if (has_ex) par[5] else NA
  u <- seq(-6, 6, length.out = ng)
  h <- u[2] - u[1]
  pz <- dnorm(u) * h
  ll <- 0
  lagm <- ar * s * u
  for (idx in split(seq_len(nrow(d)), d$g)) {
    i <- idx[1]
    f <- pz * dens(d$y[i], inv(bb0 + bb1 * d$x[i] + s * u), ex, d$cen[i])
    c0 <- sum(f)
    if (!is.finite(c0) || c0 <= 0) return(1e10)
    ll <- ll + log(c0)
    f <- f / c0
    for (i in idx[-1]) {
      eta <- bb0 + bb1 * d$x[i] + s * u
      M <- if (post) outer(lagm, inv(eta), "+") else inv(outer(lagm, eta, "+"))
      P <- matrix(dens(d$y[i], as.vector(M), ex, d$cen[i]), ng, ng)
      f <- as.vector(crossprod(f, P)) * pz
      c0 <- sum(f)
      if (!is.finite(c0) || c0 <= 0) return(1e10)
      ll <- ll + log(c0)
      f <- f / c0
    }
  }
  -ll
}
exact_nll <- if (kind == "T") exact_T else exact_F

hyst <- function(f, seed) {
  env <- f$obj$env
  par <- f$opt$par
  r <- env$random
  bb <- f$estimates$b
  vals <- numeric(0)
  nan <- 0L
  for (i in 1:12) {
    set.seed(seed * 100 + i)
    env$last.par.best[r] <- if (i == 1) 0 else if (i == 2) bb else {
      bb + rnorm(length(r), 0, 0.25 * (i %% 5))
    }
    v <- tryCatch(suppressWarnings(f$obj$fn(par)), error = function(e) NaN)
    if (is.finite(v)) vals <- c(vals, as.numeric(v)) else nan <- nan + 1L
  }
  c(nan = nan, span = if (length(vals)) diff(range(vals)) else NA_real_)
}

one <- function(seed) {
  d <- simd(seed)
  lhs <- paste0("y", if (dist == "binom") " | trials(tr)", if (cens) " | cens(cen)")
  rhs <- if (kind == "T") " ~ x + ar(t, g, cov = TRUE)" else " ~ x + ar(t, g)"
  w <- character(0); e <- ""
  f <- withCallingHandlers(tryCatch(frm(as.formula(paste0(lhs, rhs)), family = fam, data = d),
                                    error = function(x) {
                                      e <<- substr(conditionMessage(x), 1, 90); NULL
                                    }),
                           warning = function(c) {
                             w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
                           }, message = function(c) invokeRestart("muffleMessage"))
  out <- list(design = design, seed = seed, err = gsub("[\t\r\n]", " ", e), nwarn = length(w))
  st <- c(b0, b1, log(s0),
          if (kind == "T") ar0 / sqrt(1 - ar0^2) else qlogis((ar0 + 1) / 2),
          if (has_ex) log(shape0))
  if (!is.null(f)) {
    th <- f$estimates$theta
    bd <- f$opt$par[names(f$opt$par) == "betad"]
    pf <- c(f$estimates$beta[1:2], th[1:2], if (has_ex) bd[1])
    out$ll <- as.numeric(logLik(f)); out$code <- f$opt$convergence
    out$sderr <- exp(th[1])
    out$ar <- if (kind == "T") th[2] / sqrt(1 + th[2]^2) else 2 / (1 + exp(-th[2])) - 1
    out$lshape <- if (has_ex) bd[1] else NA
    vc <- tryCatch(as.matrix(vcov(f, full = TRUE)), error = function(x) NULL)
    out$se_ls <- tryCatch(sqrt(vc["theta_1", "theta_1"]), error = function(x) NA)
    h <- hyst(f, seed)
    out$nan_starts <- h[["nan"]]; out$span <- h[["span"]]
    pc <- pf
    if (kind == "T") pc[4] <- pmin(pmax(pc[4], -8), 8)
    if (has_ex) pc[5] <- pmin(pc[5], log(1e6))
    out$x_at_frm <- -exact_nll(pc, d)
  }
  lo <- c(-Inf, -Inf, log(0.01), if (kind == "T") -8 else -12, if (has_ex) log(0.05))
  hi <- c(Inf, Inf, log(10), if (kind == "T") 8 else 12, if (has_ex) log(1e4))
  ex <- tryCatch(nlminb(st, exact_nll, d = d, lower = lo, upper = hi,
                        control = list(iter.max = 60, eval.max = 150)),
                 error = function(x) NULL)
  if (!is.null(ex)) {
    out$x_ll <- -ex$objective
    out$x_sderr <- exp(ex$par[3])
    out$x_ar <- if (kind == "T") ex$par[4] / sqrt(1 + ex$par[4]^2) else 2 / (1 + exp(-ex$par[4])) - 1
  }
  out
}
cols <- c("design", "seed", "err", "nwarn", "ll", "code", "sderr", "ar", "lshape", "se_ls",
          "nan_starts", "span", "x_at_frm", "x_ll", "x_sderr", "x_ar")
first <- !file.exists(outf)
for (sd in s1:s2) {
  o <- one(sd)
  for (cn in setdiff(cols, names(o))) o[[cn]] <- NA
  write.table(as.data.frame(o[cols]), outf, sep = "\t", row.names = FALSE,
              quote = FALSE, col.names = first, append = !first)
  first <- FALSE
}
