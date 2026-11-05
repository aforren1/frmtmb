# Short check, item 2: M1's warning on designs not used to calibrate it, against
# the exact marginal by forward recursion. Each dataset records whether the
# warning fired, and the row curvature. Lognormal and the gaussian-like families
# on an identity link have a curvature in the data's units, so a rescale tests
# whether the floor is about the approximation at all.
# Usage: Rscript review-latres-m1b.R <design> <seed_from> <seed_to> <outfile>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
options(warn = -1)
args <- commandArgs(TRUE)
design <- args[[1]]
s1 <- as.integer(args[[2]])
s2 <- as.integer(args[[3]])
outf <- args[[4]]

D0 <- list(
  f_pois05 = list(fam = quote(poisson()), dist = "pois", b0 = -.7, b1 = .4, s = .7, ar = .3),
  f_pois_ar8 = list(fam = quote(poisson()), dist = "pois", b0 = .7, b1 = .4, s = .3, ar = .8),
  f_pois_short = list(fam = quote(poisson()), dist = "pois", G = 60, K = 4, b0 = .5, b1 = .4, s = .6, ar = .5),
  f_pois_long = list(fam = quote(poisson()), dist = "pois", G = 6, K = 40, b0 = 0, b1 = .4, s = .5, ar = .7),
  f_bin2 = list(fam = quote(binomial()), dist = "binom", tr = 2, b0 = 0, b1 = .5, s = .7, ar = .5),
  f_bin6 = list(fam = quote(binomial()), dist = "binom", tr = 6, b0 = -1, b1 = .5, s = .6, ar = .5),
  f_bin30 = list(fam = quote(binomial()), dist = "binom", tr = 30, b0 = -1.5, b1 = .4, s = .5, ar = .6),
  f_probit = list(fam = quote(bernoulli("probit")), dist = "bern", b0 = 0, b1 = .4, s = .6, ar = .6),
  f_cloglog = list(fam = quote(bernoulli("cloglog")), dist = "bern", b0 = -.5, b1 = .4, s = .8, ar = .5),
  f_geom = list(fam = quote(geometric()), dist = "geom", b0 = 1, b1 = .4, s = .5, ar = .6),
  f_nb1 = list(fam = quote(negbinomial()), dist = "nb", shape = 1, b0 = 1, b1 = .4, s = .5, ar = .6),
  f_nb8 = list(fam = quote(negbinomial()), dist = "nb", shape = 8, b0 = 2, b1 = .4, s = .8, ar = .4),
  f_gam1 = list(fam = quote(Gamma(link = "log")), dist = "gamma", shape = 1, b0 = 0, b1 = .4, s = .5, ar = .6),
  f_gam2s1 = list(fam = quote(Gamma(link = "log")), dist = "gamma", shape = 2, b0 = 0, b1 = .4, s = 1, ar = .5),
  f_gam20 = list(fam = quote(Gamma(link = "log")), dist = "gamma", shape = 20, b0 = 0, b1 = .4, s = .3, ar = .7),
  f_wei2 = list(fam = quote(weibull()), dist = "wei", shape = 2, b0 = 0, b1 = .4, s = .5, ar = .5),
  f_exp = list(fam = quote(exponential()), dist = "expo", b0 = 0, b1 = .4, s = .7, ar = .5),
  f_beta2 = list(fam = quote(Beta()), dist = "beta", phi = 2, b0 = 0, b1 = .4, s = .5, ar = .5),
  f_beta20 = list(fam = quote(Beta()), dist = "beta", phi = 20, b0 = 0, b1 = .4, s = .4, ar = .7),
  f_zinb = list(fam = quote(zero_inflated_negbinomial()), dist = "zinb", zi = .2, shape = 3, b0 = 1.2, b1 = .4, s = .5, ar = .5),
  f_hnb = list(fam = quote(hurdle_negbinomial()), dist = "hnb", hu = .3, shape = 3, b0 = 1, b1 = .4, s = .5, ar = .5),
  f_cumpro = list(fam = quote(cumulative("probit")), dist = "cum", pfun = pnorm, rfun = rnorm, thr = c(-1, 0, 1), b0 = 0, b1 = .4, s = .6, ar = .5),
  f_ln1 = list(fam = quote(lognormal()), dist = "ln", sigma = 1, b0 = 0, b1 = .4, s = .5, ar = .5))
D <- D0[[design]]
fam <- eval(D$fam, asNamespace("frmtmb"))
linkinv <- fam$linkinv %||% fam$links$mu$linkinv
G <- D$G %||% 25
K <- D$K %||% 8
dsn <- function(y, xi, omega, alpha) 2 / omega * dnorm((y - xi) / omega) * pnorm(alpha * (y - xi) / omega)
sn_xi <- function(mu, omega, alpha) mu - omega * (alpha / sqrt(1 + alpha^2)) * sqrt(2 / pi)
dens <- function(y, eta, p) {
  mu <- linkinv(eta)
  v <- switch(D$dist,
    pois = dpois(y, mu), bern = dbinom(y, 1, mu), geom = dnbinom(y, size = 1, mu = mu),
    gamma = dgamma(y, p$shape, p$shape / mu), wei = dweibull(y, p$shape, mu / gamma(1 + 1 / p$shape)), expo = dexp(y, 1 / mu),
    zinb = if (y == 0) p$zi + (1 - p$zi) * dnbinom(0, size = p$shape, mu = mu) else (1 - p$zi) * dnbinom(y, size = p$shape, mu = mu),
    hnb = if (y == 0) p$hu + 0 * mu else (1 - p$hu) * dnbinom(y, size = p$shape, mu = mu) / (1 - dnbinom(0, size = p$shape, mu = mu)),
    zip = if (y == 0) p$zi + (1 - p$zi) * dpois(0, mu) else (1 - p$zi) * dpois(y, mu),
    hp = if (y == 0) p$hu + 0 * mu else (1 - p$hu) * dpois(y, mu) / (1 - dpois(0, mu)),
    zib = if (y == 0) p$zi + (1 - p$zi) * dbinom(0, D$tr, mu) else (1 - p$zi) * dbinom(y, D$tr, mu),
    cum = {
      th <- c(-Inf, p$thr, Inf)
      (D$pfun %||% plogis)(th[y + 1] - eta) - (D$pfun %||% plogis)(th[y] - eta)
    },
    nb = dnbinom(y, size = p$shape, mu = mu),
    beta = dbeta(y, mu * p$phi, (1 - mu) * p$phi),
    binom = dbinom(y, D$tr, mu),
    ln = dlnorm(y, eta, p$sigma),
    sn = dsn(y, sn_xi(eta, p$sigma, p$alpha), p$sigma, p$alpha))
  v[!is.finite(v)] <- 0
  v
}
simd <- function(seed) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = seq_len(G))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- rnorm(n)
  e <- as.vector(replicate(G, arima.sim(list(ar = D$ar), K, sd = D$s)))
  eta <- D$b0 + D$b1 * d$x + e
  mu <- linkinv(eta)
  d$y <- switch(D$dist,
    pois = rpois(n, mu), bern = rbinom(n, 1, mu), geom = rnbinom(n, size = 1, mu = mu),
    gamma = rgamma(n, D$shape, D$shape / mu), wei = rweibull(n, D$shape, mu / gamma(1 + 1 / D$shape)), expo = rexp(n, 1 / mu),
    zinb = rnbinom(n, size = D$shape, mu = mu) * rbinom(n, 1, 1 - D$zi),
    hnb = { pos <- rbinom(n, 1, 1 - D$hu); pmax(qnbinom(runif(n, dnbinom(0, size = D$shape, mu = mu), 1), size = D$shape, mu = mu), 1) * pos },
    zip = rpois(n, mu) * rbinom(n, 1, 1 - D$zi),
    hp = {
      pos <- rbinom(n, 1, 1 - D$hu)
      yy <- qpois(runif(n, dpois(0, mu), 1), mu)
      pmax(yy, 1) * pos
    },
    zib = rbinom(n, D$tr, mu) * rbinom(n, 1, 1 - D$zi),
    cum = findInterval(eta + (D$rfun %||% rlogis)(n), D$thr) + 1L,
    nb = rnbinom(n, size = D$shape, mu = mu),
    beta = pmin(pmax(rbeta(n, mu * D$phi, (1 - mu) * D$phi), 1e-6), 1 - 1e-6),
    binom = rbinom(n, D$tr, mu),
    ln = rlnorm(n, eta, D$sigma),
    sn = {
      dl <- D$alpha / sqrt(1 + D$alpha^2)
      sn_xi(mu, D$sigma, D$alpha) + D$sigma * (dl * abs(rnorm(n)) + sqrt(1 - dl^2) * rnorm(n))
    })
  if (D$dist == "cum") d$y <- factor(d$y, ordered = TRUE)
  d$tr <- D$tr %||% 1L
  d$g <- factor(d$g)
  d
}
extra_par <- function(p) {
  switch(D$dist,
    pois = , bern = , geom = , expo = list(), gamma = , wei = list(shape = exp(p[1])),
    zinb = list(zi = plogis(p[1]), shape = exp(p[2])), hnb = list(hu = plogis(p[1]), shape = exp(p[2])),
    zip = , zib = list(zi = plogis(p[1])), hp = list(hu = plogis(p[1])),
    cum = list(thr = c(p[1], p[1] + cumsum(exp(p[-1])))),
    nb = list(shape = exp(p[1])), beta = list(phi = exp(p[1])),
    binom = list(), ln = list(sigma = exp(p[1])),
    sn = list(sigma = exp(p[1]), alpha = D$alpha))
}
extra_start <- function() {
  switch(D$dist,
    pois = , bern = , geom = , expo = numeric(0), gamma = , wei = log(D$shape),
    zinb = c(qlogis(D$zi), log(D$shape)), hnb = c(qlogis(D$hu), log(D$shape)),
    zip = , zib = qlogis(D$zi), hp = qlogis(D$hu),
    cum = c(D$thr[1], log(diff(D$thr))),
    nb = log(D$shape), beta = log(D$phi), binom = numeric(0), ln = log(D$sigma), sn = log(D$sigma))
}
is_cum <- D$dist == "cum"
exact <- function(par, d, ng = 121) {
  # cumulative: no intercept in the fit, the thresholds carry it
  b0 <- if (is_cum) 0 else par[1]
  b1 <- par[2]
  s <- exp(par[3])
  th <- par[4]
  p <- extra_par(par[-(1:4)])
  msd <- s * sqrt(1 + th^2)
  ar <- th / sqrt(1 + th^2)
  u <- seq(-6.5, 6.5, length.out = ng)
  h <- u[2] - u[1]
  grid <- u * msd
  Tm <- outer(grid, grid, function(a, b) dnorm(b, ar * a, s)) * (h * msd)
  y <- if (is_cum) as.integer(d$y) else d$y
  ll <- 0
  for (idx in split(seq_len(nrow(d)), d$g)) {
    f <- dnorm(u) * h
    for (k in seq_along(idx)) {
      i <- idx[k]
      if (k > 1L) f <- as.vector(crossprod(f, Tm))
      f <- f * dens(y[i], b0 + b1 * d$x[i] + grid, p)
      c0 <- sum(f)
      if (!is.finite(c0) || c0 <= 0) return(1e10)
      ll <- ll + log(c0)
      f <- f / c0
    }
  }
  -ll
}
one <- function(seed) {
  d <- simd(seed)
  lhs <- if (D$dist == "binom" || D$dist == "zib") "y | trials(tr)" else "y"
  form <- as.formula(paste(lhs, "~ x + ar(t, g, cov = TRUE)"))
  w <- character(0)
  f <- withCallingHandlers(
    tryCatch(frm(form, data = d, family = fam), error = function(e) e),
    warning = function(c) {
      w <<- c(w, paste(class(c)[1], conditionMessage(c)))
      invokeRestart("muffleWarning")
    }, message = function(c) invokeRestart("muffleMessage"))
  out <- list(design = design, seed = seed)
  if (inherits(f, "error")) {
    out$err <- substr(conditionMessage(f), 1, 90)
  } else {
    out$err <- ""
    th <- f$estimates$theta
    out$ll <- as.numeric(logLik(f)); out$code <- f$opt$convergence
    out$b1 <- f$estimates$beta[if (is_cum) 1 else 2]
    out$ls <- th[1]
    out$t <- th[2]
    out$ar <- th[2] / sqrt(1 + th[2]^2)
    out$warned <- as.integer(any(grepl("Laplace approximation of", w)))
    out$info <- stats::median(fi("latres_row_info")(f, fi("frame_latres")(f$frame)[[1]]))
    vc <- tryCatch(as.matrix(vcov(f, full = TRUE)), error = function(e) NULL)
    out$se_ls <- tryCatch(sqrt(vc["theta_1", "theta_1"]), error = function(e) NA)
    out$se_b1 <- tryCatch(sqrt(vc[if (is_cum) "x" else 2, if (is_cum) "x" else 2]), error = function(e) NA)
  }
  st <- c(D$b0, D$b1, log(D$s), D$ar / sqrt(1 - D$ar^2), extra_start())
  lo <- c(-Inf, -Inf, log(.01), -8, rep(-Inf, length(st) - 4))
  hi <- c(Inf, Inf, log(30), 8, rep(Inf, length(st) - 4))
  fnx <- if (is_cum) function(p, d) exact(c(0, p[-1]), d) else exact
  ex <- tryCatch(nlminb(st, fnx, d = d, lower = lo, upper = hi, control = list(iter.max = 80, eval.max = 200)), error = function(e) NULL)
  if (!is.null(ex)) {
    out$x_b1 <- ex$par[2]; out$x_ls <- ex$par[3]; out$x_t <- ex$par[4]
    out$x_ar <- ex$par[4] / sqrt(1 + ex$par[4]^2)
    out$x_ll <- -ex$objective
  }
  out
}
cols <- c("design","seed","err","ll","code","b1","ls","t","ar","warned","info","se_ls","se_b1","x_b1","x_ls","x_t","x_ar","x_ll")
first <- TRUE
for (sd in s1:s2) {
  o <- one(sd)
  for (cn in setdiff(cols, names(o))) o[[cn]] <- NA
  write.table(as.data.frame(o[cols]), outf, sep = "\t", row.names = FALSE, quote = TRUE, col.names = first, append = !first)
  first <- FALSE
}
