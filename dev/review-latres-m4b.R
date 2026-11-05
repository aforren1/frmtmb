# Short check, item 1: is the refusal of M4 the right size? Benign designs, a
# mean far from its bound and a low sderr, read past the refusal (in-process
# no-op), against the exact marginal by forward recursion. For each fit: the
# distance of the mean from its bound, in sds of the mean, as a candidate for a
# condition on the fit.
# Usage: Rscript review-latres-m4b.R <design> <seed_from> <seed_to> <outfile>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
args <- commandArgs(TRUE)
design <- args[[1]]
s1 <- as.integer(args[[2]])
s2 <- as.integer(args[[3]])
outf <- args[[4]]
assignInNamespace("latres_refuse_laplace", function(spec, frame) invisible(NULL), ns = "frmtmb"); options(warn = -1)

dig <- function(y, mu, shape) sqrt(shape / (2 * pi * y^3)) * exp(-shape * (y - mu)^2 / (2 * mu^2 * y))
rig <- function(n, mu, shape) {
  v <- rnorm(n)^2
  x <- mu + mu^2 * v / (2 * shape) - mu / (2 * shape) * sqrt(4 * mu * shape * v + mu^2 * v^2)
  ifelse(runif(n) <= mu / (mu + x), x, mu^2 / x)
}
D0 <- list(
  # link scale, a link that does not hold the mean (cov = TRUE)
  ginv_b1 = list(kind = "T", fam = quote(Gamma()), dist = "gamma", shape = 5, b0 = 2, b1 = .1, s = .1, ar = .5),
  ginv_b2 = list(kind = "T", fam = quote(Gamma()), dist = "gamma", shape = 20, b0 = 2, b1 = .1, s = .2, ar = .5),
  ginv_b3 = list(kind = "T", fam = quote(Gamma()), dist = "gamma", shape = 50, b0 = 2, b1 = .1, s = .3, ar = .5),
  gident_b = list(kind = "T", fam = quote(Gamma(link = "identity")), dist = "gamma", shape = 20, b0 = 20, b1 = 1, s = 2, ar = .5),
  ig_b1 = list(kind = "T", fam = quote(inverse.gaussian()), dist = "ig", shape = 20, b0 = 4, b1 = .2, s = .3, ar = .5),
  ig_b2 = list(kind = "T", fam = quote(inverse.gaussian()), dist = "ig", shape = 50, b0 = 4, b1 = .2, s = .5, ar = .5),
  pid_b = list(kind = "T", fam = quote(poisson("identity")), dist = "pois", b0 = 50, b1 = 3, s = 4, ar = .5),
  pid_8 = list(kind = "T", fam = quote(poisson("identity")), dist = "pois", b0 = 8, b1 = .5, s = 1.2, ar = .5),
  bid_b = list(kind = "T", fam = quote(binomial("identity")), dist = "binom", tr = 50, b0 = .5, b1 = .02, s = .04, ar = .5),
  # the same designs on a link that holds the mean: the control
  glog_c = list(kind = "T", fam = quote(Gamma(link = "log")), dist = "gamma", shape = 20, b0 = .7, b1 = .05, s = .1, ar = .5),
  # post-link, cov = FALSE, a mean far from the bound
  gpost_b = list(kind = "F", post = TRUE, inv = exp, fam = quote(Gamma(link = "log")), dist = "gamma", shape = 20, b0 = 3, b1 = .2, s = .5, ar = .5),
  wpost_b = list(kind = "F", post = TRUE, inv = exp, fam = quote(weibull()), dist = "weibull", shape = 4, b0 = 3, b1 = .2, s = .4, ar = .5),
  epost_b = list(kind = "F", post = TRUE, inv = exp, fam = quote(exponential()), dist = "exp", b0 = 3, b1 = .2, s = .4, ar = .5),
  igpost_b = list(kind = "F", post = TRUE, inv = exp, fam = quote(inverse.gaussian(link = "log")), dist = "ig", shape = 200, b0 = 3, b1 = .2, s = .4, ar = .5),
  bepost_b = list(kind = "F", post = TRUE, inv = plogis, fam = quote(Beta()), dist = "beta", phi = 50, b0 = 0, b1 = .1, s = .1, ar = .5),
  psqrt_b = list(kind = "F", post = TRUE, inv = function(e) e^2, fam = quote(poisson("sqrt")), dist = "pois", b0 = 7, b1 = .3, s = .5, ar = .5),
  # the two families the lane did not measure
  bb_post_lo = list(kind = "F", post = TRUE, inv = plogis, fam = quote(beta_binomial()), dist = "bb", tr = 10, phi = 10, b0 = 0, b1 = .2, s = .5, ar = .5),
  bb_post_b = list(kind = "F", post = TRUE, inv = plogis, fam = quote(beta_binomial()), dist = "bb", tr = 50, phi = 50, b0 = 0, b1 = .1, s = .1, ar = .5),
  zoib_post_lo = list(kind = "F", post = TRUE, inv = plogis, fam = quote(zero_one_inflated_beta()), dist = "zoib", phi = 10, zoi = .1, coi = .5, b0 = 0, b1 = .2, s = .5, ar = .5),
  zoib_post_b = list(kind = "F", post = TRUE, inv = plogis, fam = quote(zero_one_inflated_beta()), dist = "zoib", phi = 50, zoi = .1, coi = .5, b0 = 0, b1 = .1, s = .1, ar = .5),
  # lane's own harmful designs, for the margin; and weak-signal controls on a held link
  gid_l = list(kind = "T", fam = quote(Gamma(link = "identity")), dist = "gamma", shape = 4, b0 = 8, b1 = 1, s = 1, ar = .5),
  bid_l = list(kind = "T", fam = quote(binomial("identity")), dist = "binom", tr = 10, b0 = .5, b1 = .05, s = .05, ar = .5),
  ginv_l = list(kind = "T", fam = quote(Gamma()), dist = "gamma", shape = 4, b0 = 1.5, b1 = .1, s = .15, ar = .5),
  pid_l = list(kind = "T", fam = quote(poisson("identity")), dist = "pois", b0 = 2, b1 = .5, s = 1, ar = .5),
  ginv_s5 = list(kind = "T", fam = quote(Gamma()), dist = "gamma", shape = 5, b0 = 2, b1 = .1, s = .3, ar = .5),
  gid_s5 = list(kind = "T", fam = quote(Gamma(link = "identity")), dist = "gamma", shape = 5, b0 = 20, b1 = 1, s = 2, ar = .5),
  ig_s5 = list(kind = "T", fam = quote(inverse.gaussian()), dist = "ig", shape = 5, b0 = 4, b1 = .2, s = .5, ar = .5),
  glog_s5 = list(kind = "T", fam = quote(Gamma(link = "log")), dist = "gamma", shape = 5, b0 = .7, b1 = .05, s = .15, ar = .5),
  pid_20 = list(kind = "T", fam = quote(poisson("identity")), dist = "pois", b0 = 20, b1 = 1, s = 2, ar = .5),
  glog_w = list(kind = "T", fam = quote(Gamma(link = "log")), dist = "gamma", shape = 5, b0 = .7, b1 = .05, s = .05, ar = .5),
  blog_w = list(kind = "T", fam = quote(binomial()), dist = "binom", tr = 10, b0 = 0, b1 = .2, s = .2, ar = .5),
  # and held forms of those two under cov = TRUE: the lane keeps them
  bb_T = list(kind = "T", fam = quote(beta_binomial()), dist = "bb", tr = 10, phi = 10, b0 = 0, b1 = .2, s = .5, ar = .5),
  zoib_T = list(kind = "T", fam = quote(zero_one_inflated_beta()), dist = "zoib", phi = 10, zoi = .1, coi = .5, b0 = 0, b1 = .2, s = .5, ar = .5))
D <- D0[[design]]
fam <- eval(D$fam, asNamespace("frmtmb"))
G <- 25
K <- 6
ngrid <- if (D$kind == "F") 41 else 121
linkinv <- fam$linkinv %||% fam$links$mu$linkinv

dens <- function(y, mu, tr = D$tr) {
  shape <- D$shape; phi <- D$phi
  v <- switch(D$dist,
    gamma = dgamma(y, shape, shape / mu),
    ig = dig(y, mu, shape),
    pois = dpois(y, mu),
    binom = dbinom(y, tr, mu),
    weibull = dweibull(y, shape, mu / gamma(1 + 1 / shape)),
    exp = dexp(y, 1 / mu),
    beta = dbeta(y, mu * phi, (1 - mu) * phi),
    bb = exp(lchoose(tr, y) + lbeta(y + mu * phi, tr - y + (1 - mu) * phi) - lbeta(mu * phi, (1 - mu) * phi)),
    zoib = if (y == 0) D$zoi * (1 - D$coi) + 0 * mu else if (y == 1) D$zoi * D$coi + 0 * mu else (1 - D$zoi) * dbeta(y, mu * phi, (1 - mu) * phi))
  bad <- switch(D$dist, gamma = , ig = , pois = , weibull = , exp = mu <= 0, binom = , beta = , bb = , zoib = mu <= 0 | mu >= 1)
  v[bad | !is.finite(v)] <- 0
  v
}
# the design's own dispersion: the exact search moves D's copy
D_true <- D
draw <- function(mu) {
  n <- length(mu)
  shape <- D_true$shape; phi <- D_true$phi
  switch(D$dist,
    gamma = rgamma(n, shape, shape / mu), ig = rig(n, mu, shape), pois = rpois(n, mu),
    binom = rbinom(n, D$tr, mu), weibull = rweibull(n, shape, mu / gamma(1 + 1 / shape)),
    exp = rexp(n, 1 / mu), beta = rbeta(n, mu * phi, (1 - mu) * phi),
    bb = rbinom(n, D$tr, rbeta(n, mu * phi, (1 - mu) * phi)),
    zoib = {
      u <- runif(n)
      ifelse(u < D_true$zoi, rbinom(n, 1, D_true$coi), rbeta(n, mu * phi, (1 - mu) * phi))
    })
}
okmu <- function(mu) switch(D$dist, binom = , beta = , bb = , zoib = all(mu > 0 & mu < 1), all(mu > 0))
simd <- function(seed) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = seq_len(G))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  d$y <- NA_real_
  for (g in seq_len(G)) {
    i <- which(d$g == g)
    repeat {
      if (D$kind == "T") {
        mu <- linkinv(D$b0 + D$b1 * d$x[i] + as.numeric(arima.sim(list(ar = D$ar), K, sd = D$s)))
      } else {
        z <- rnorm(K)
        mu <- D$inv(D$b0 + D$b1 * d$x[i] + D$s * z) + D$ar * D$s * c(0, z[-K])
      }
      if (okmu(mu)) break
    }
    d$y[i] <- draw(mu)
  }
  d$tr <- D$tr %||% 1L
  d$g <- factor(d$g)
  d
}
# exact: b0, b1, log s, chart of ar; dispersion as frm() holds it is fixed at the
# truth in the exact search where it would add a dimension to nothing: it is free
npar <- 4 + (!is.null(D$shape)) + (!is.null(D$phi))
set_disp <- function(p) {
  if (!is.null(D$shape)) D$shape <<- exp(p[5])
  if (!is.null(D$phi)) D$phi <<- exp(p[5])
}
exact_T <- function(par, d, ng = ngrid) {
  b0 <- par[1]; b1 <- par[2]; s <- exp(par[3]); th <- par[4]
  if (length(par) > 4) set_disp(par)
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
      f <- f * dens(d$y[i], linkinv(b0 + b1 * d$x[i] + grid))
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
  if (length(par) > 4) set_disp(par)
  u <- seq(-6, 6, length.out = ng)
  pz <- dnorm(u) * (u[2] - u[1])
  lagm <- ar * s * u
  ll <- 0
  for (idx in split(seq_len(nrow(d)), d$g)) {
    f <- pz * dens(d$y[idx[1]], D$inv(b0 + b1 * d$x[idx[1]] + s * u))
    c0 <- sum(f)
    if (!is.finite(c0) || c0 <= 0) return(1e10)
    ll <- ll + log(c0)
    f <- f / c0
    for (i in idx[-1]) {
      M <- outer(lagm, D$inv(b0 + b1 * d$x[i] + s * u), "+")
      P <- matrix(dens(d$y[i], as.vector(M)), ng, ng)
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
truth_disp <- D$shape %||% D$phi

one <- function(seed) {
  d <- simd(seed)
  lhs <- if (D$dist %in% c("binom", "bb")) "y | trials(tr)" else "y"
  form <- as.formula(paste(lhs, if (D$kind == "T") "~ x + ar(t, g, cov = TRUE)" else "~ x + ar(t, g)"))
  m <- msg_of(frm(form, data = d, family = fam))
  out <- list(design = design, seed = seed, err = substr(m$err %||% "", 1, 80))
  dp <- NULL
  if (is.null(m$err)) {
    f <- m$val
    par <- f$opt$par
    nm <- names(par)
    th <- f$estimates$theta
    ib <- which(nm == "betad")
    dp <- if (length(ib)) par[ib[1]] else NA
    out$ll <- as.numeric(logLik(f)); out$code <- f$opt$convergence
    out$b1 <- f$estimates$beta[2]; out$ls <- th[1]; out$t <- th[2]
    out$ar <- if (D$kind == "T") th[2] / sqrt(1 + th[2]^2) else 2 / (1 + exp(-th[2])) - 1
    out$ldisp <- dp
    vc <- tryCatch(as.matrix(vcov(f, full = TRUE)), error = function(e) NULL)
    out$se_b1 <- tryCatch(sqrt(vc[2, 2]), error = function(e) NA)
    out$se_ls <- tryCatch(sqrt(vc["theta_1", "theta_1"]), error = function(e) NA)
    out$se_t <- tryCatch(sqrt(vc["theta_2", "theta_2"]), error = function(e) NA)
    out$nwarn <- length(m$warn)
    # the candidate condition: how far the mean stays from its bound, in sds
    # of the mean's latent part, at the fit
    bk <- fi("frame_latres")(f$frame)[[1]]
    fitv <- as.numeric(as.matrix(fitted(f))[, 1])
    s <- exp(th[1])
    a <- bk$aux
    msd <- fi("latres_marginal_sd")(th, a)
    et <- as.numeric(frm_linpred(f, type = "link"))
    dist_b <- switch(D$dist, bb = pmin(fitv / D$tr, 1 - fitv / D$tr), binom = , beta = , zoib = pmin(fitv, 1 - fitv), fitv)
    lnk <- fam$link %||% fam$links$mu$name %||% "identity"
    if (D$kind == "T") {
      # the way out of the support is eta crossing the link's own edge
      if (lnk %in% c("identity", "inverse", "1/mu^2")) {
        de <- if (D$dist %in% c("binom", "beta", "bb", "zoib")) pmin(et, 1 - et) else et
        out$margin_min <- min(de / msd)
        out$margin_med <- median(de / msd)
      }
    } else {
      # the way out is the AR sum after the inverse link
      out$margin_min <- min(dist_b / (abs(out$ar) * s))
      out$margin_med <- median(dist_b / (abs(out$ar) * s))
    }
    if (D$dist == "zoib") {
      bd <- par[names(par) == "betad"]
      D$zoi <<- plogis(bd[2])
      D$coi <<- plogis(bd[3])
    }
    # the objective over inner starts
    env <- f$obj$env
    r <- env$random
    bb <- f$estimates$b
    vals <- numeric(0)
    nan <- 0L
    for (i in 1:8) {
      set.seed(seed * 100 + i)
      env$last.par.best[r] <- if (i == 1) 0 else bb + rnorm(length(r), 0, .1 * i)
      v <- tryCatch(suppressWarnings(as.numeric(f$obj$fn(par))), error = function(e) NaN)
      if (is.finite(v)) vals <- c(vals, v) else nan <- nan + 1L
    }
    out$span <- if (length(vals)) max(vals) - min(vals) else NA
    out$nan <- nan
    pf <- c(out$b0 <- f$estimates$beta[1], out$b1, th[1], if (D$kind == "T") pmin(pmax(th[2], -8), 8) else pmin(pmax(th[2], -12), 12), if (npar > 4) min(dp, log(1e5)))
    out$x_at_frm <- -exact_nll(pf, d)
  }
  st <- c(D$b0, D$b1, log(D$s), if (D$kind == "T") D$ar / sqrt(1 - D$ar^2) else qlogis((D$ar + 1) / 2), if (npar > 4) log(truth_disp))
  lo <- c(-Inf, -Inf, log(0.003), if (D$kind == "T") -8 else -12, if (npar > 4) log(.05))
  hi <- c(Inf, Inf, log(20), if (D$kind == "T") 8 else 12, if (npar > 4) log(1e5))
  ex <- tryCatch(nlminb(st, exact_nll, d = d, lower = lo, upper = hi, control = list(iter.max = 60, eval.max = 150)), error = function(e) NULL)
  if (!is.null(ex)) {
    out$x_b1 <- ex$par[2]; out$x_ls <- ex$par[3]; out$x_t <- ex$par[4]
    out$x_ar <- if (D$kind == "T") ex$par[4] / sqrt(1 + ex$par[4]^2) else 2 / (1 + exp(-ex$par[4])) - 1
    out$x_ldisp <- if (npar > 4) ex$par[5] else NA
    out$x_ll <- -ex$objective
  }
  out
}
cols <- c("design","seed","err","ll","code","b0","b1","ls","t","ar","ldisp","se_b1","se_ls","se_t","nwarn","margin_min","margin_med","span","nan","x_at_frm","x_b1","x_ls","x_t","x_ar","x_ldisp","x_ll")
first <- TRUE
for (sd in s1:s2) {
  o <- one(sd)
  for (cn in setdiff(cols, names(o))) o[[cn]] <- NA
  write.table(as.data.frame(o[cols]), outf, sep = "\t", row.names = FALSE, quote = TRUE, col.names = first, append = !first)
  first <- FALSE
}
