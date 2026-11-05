# M1: which designs does Laplace read wrongly under a latent residual, and what
# statistic of the fit separates them? ar(t, g, cov = TRUE), the exact MLE by the
# forward recursion of dev/review-latres-sim.R, adapted (families added; the
# same chain and grid rule, 121 points over +-6.5 marginal sds).
#
# Usage: Rscript dev/latres-m1.R <lib> <design> <seed_from> <seed_to> <outfile>
# Per dataset: frm(), the statistic candidates, the exact MLE. Data seed = seed.
args <- commandArgs(TRUE)
.libPaths(c(args[[1]], "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
design <- args[[2]]; s1 <- as.integer(args[[3]]); s2 <- as.integer(args[[4]])
outf <- args[[5]]
`%||%` <- function(a, b) if (is.null(a)) b else a
P <- function(G, K, b0, b1, s, ar) list(G = G, K = K, b0 = b0, b1 = b1, s = s, ar = ar)
dz <- list(
  # the reviewer's
  bern = c(fam = "bernoulli", P(20, 10, 0.2, 0.5, 1, 0.6)),
  bern_long = c(fam = "bernoulli", P(5, 60, 0.2, 0.5, 1, 0.8)),
  binom3 = c(fam = "binomial", tr = 3, P(15, 10, 0.2, 0.5, 0.8, 0.6)),
  pois = c(fam = "poisson", P(10, 15, -0.2, 0.5, 0.5, 0.6)),
  pois_ar9 = c(fam = "poisson", P(10, 15, -0.2, 0.5, 0.3, 0.9)),
  pois_s1 = c(fam = "poisson", P(10, 15, -0.5, 0.5, 1, 0.6)),
  pois_big = c(fam = "poisson", P(10, 15, 3, 0.3, 0.5, 0.6)),
  # mine: the trials ladder, weak latents, low counts, longer series
  bern_s05 = c(fam = "bernoulli", P(20, 10, 0.2, 0.5, 0.5, 0.6)),
  bern_big = c(fam = "bernoulli", P(60, 10, 0.2, 0.5, 1, 0.6)),
  binom5 = c(fam = "binomial", tr = 5, P(15, 10, 0.2, 0.5, 0.8, 0.6)),
  binom10 = c(fam = "binomial", tr = 10, P(15, 10, 0.2, 0.5, 0.8, 0.6)),
  binom20 = c(fam = "binomial", tr = 20, P(15, 10, 0.2, 0.5, 0.8, 0.6)),
  binom10_lo = c(fam = "binomial", tr = 10, P(15, 10, -2.5, 0.5, 0.8, 0.6)),
  pois_lo = c(fam = "poisson", P(10, 15, -1.5, 0.5, 0.5, 0.6)),
  pois_vlo = c(fam = "poisson", P(20, 15, -2.5, 0.5, 0.8, 0.6)),
  pois_3 = c(fam = "poisson", P(10, 15, 1, 0.5, 0.5, 0.6)),
  nb = c(fam = "negbinomial", shape = 4, P(15, 10, 1.5, 0.4, 0.5, 0.5)),
  nb_lo = c(fam = "negbinomial", shape = 4, P(15, 10, -0.2, 0.4, 0.5, 0.5)),
  gamma = c(fam = "gamma", shape = 4, P(15, 10, 0.5, 0.4, 0.5, 0.5)),
  ord4 = c(fam = "cumulative", P(20, 10, 0, 0.5, 1, 0.6)))
D <- dz[[design]]
famf <- switch(D$fam,
  bernoulli = bernoulli(), binomial = binomial(), poisson = poisson(),
  negbinomial = negbinomial(), gamma = Gamma(link = "log"),
  cumulative = cumulative())
cuts <- c(-1, 0.5, 2)

simd <- function(seed) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(D$K), g = seq_len(D$G))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  e <- numeric(nrow(d))
  for (g in seq_len(D$G)) {
    idx <- which(d$g == g)
    e[idx[1]] <- rnorm(1, 0, D$s / sqrt(1 - D$ar^2))
    for (i in idx[-1]) e[i] <- D$ar * e[i - 1] + rnorm(1, 0, D$s)
  }
  eta <- D$b0 + D$b1 * d$x + e
  n <- nrow(d)
  d$y <- switch(D$fam,
    bernoulli = rbinom(n, 1, plogis(eta)),
    binomial = rbinom(n, D$tr, plogis(eta)),
    poisson = rpois(n, exp(eta)),
    negbinomial = rnbinom(n, size = D$shape, mu = exp(eta)),
    gamma = rgamma(n, D$shape, D$shape / exp(eta)),
    cumulative = findInterval(eta + rlogis(n), cuts) + 1L)
  if (D$fam == "cumulative") d$y <- factor(d$y, ordered = TRUE)
  if (D$fam == "binomial") d$tr <- D$tr
  d$g <- factor(d$g)
  d
}
row_dens <- function(y, eta, ex) {
  switch(D$fam,
    bernoulli = dbinom(y, 1, plogis(eta)),
    binomial = dbinom(y, D$tr, plogis(eta)),
    poisson = dpois(y, exp(eta)),
    negbinomial = dnbinom(y, size = exp(ex), mu = exp(eta)),
    gamma = dgamma(y, exp(ex), exp(ex) / exp(eta)),
    cumulative = {
      th <- c(-Inf, ex, Inf)
      plogis(th[y + 1L] - eta) - plogis(th[y] - eta)
    })
}
# par: b1, s chart, ar chart, then b0 or the thresholds, then a log dispersion
n_extra <- switch(D$fam, cumulative = 3L, negbinomial = , gamma = 2L, 1L)
exact_nll <- function(par, d, ngrid = 121) {
  b1 <- par[1]; s <- exp(par[2]); th <- par[3]
  ar <- th / sqrt(1 + th^2)
  msd <- s * sqrt(1 + th^2)
  rest <- par[-(1:3)]
  b0 <- if (D$fam == "cumulative") 0 else rest[1]
  ex <- if (D$fam == "cumulative") rest else rest[2]
  if (D$fam == "cumulative" && is.unsorted(ex)) return(1e10)
  u <- seq(-6.5, 6.5, length.out = ngrid)
  h <- u[2] - u[1]
  grid <- u * msd
  Tm <- outer(grid, grid, function(a, b) dnorm(b, ar * a, s)) * (h * msd)
  ll <- 0
  yy <- if (is.factor(d$y)) as.integer(d$y) else d$y
  for (idx in split(seq_len(nrow(d)), d$g)) {
    f <- dnorm(u) * h
    for (k in seq_along(idx)) {
      i <- idx[k]
      if (k > 1L) f <- as.vector(crossprod(f, Tm))
      f <- f * row_dens(yy[i], b0 + b1 * d$x[i] + grid, ex)
      c0 <- sum(f)
      if (!is.finite(c0) || c0 <= 0) return(1e10)
      ll <- ll + log(c0)
      f <- f / c0
    }
  }
  -ll
}
# curvature of each row's log density in its linear predictor, by a second
# difference at the fit's own mean: the information one row gives its latent
row_info <- function(fit) {
  bk <- ns$frame_latres(fit$frame)[[1L]]
  r <- bk$aux$ac$resp
  rspec <- fit$spec$responses[[r]]
  fam <- rspec$family
  dp <- ns$eval_dpars(fit)[[r]]
  lk <- fit$frame$linpreds[[paste0(r, ".mu")]][["link"]] %||%
    fit$frame$linpreds[[1L]][["link"]]
  eta <- lk$linkfun(as.numeric(dp$mu))
  y <- fit$frame$y[[r]]
  av <- fit$frame$aterm_values[[r]]
  ex <- ns$fit_extras(fit)
  ll <- function(e) {
    dpe <- dp
    dpe$mu <- lk$linkinv(e)
    as.numeric(ns$row_lpdf(fam, y, y, dpe, av, ex))
  }
  h <- 1e-3
  pmax(-(ll(eta + h) - 2 * ll(eta) + ll(eta - h)) / h^2, 0)
}
one <- function(seed) {
  d <- simd(seed)
  lhs <- if (D$fam == "binomial") "y | trials(tr)" else "y"
  form <- as.formula(paste(lhs, "~ x + ar(t, g, cov = TRUE)"))
  w <- character(0)
  fit <- withCallingHandlers(tryCatch(frm(form, data = d, family = famf),
                                      error = function(e) NULL),
                             warning = function(c) {
                               w <<- c(w, conditionMessage(c))
                               invokeRestart("muffleWarning")
                             }, message = function(c) invokeRestart("muffleMessage"))
  out <- list(design = design, seed = seed)
  if (is.null(fit)) return(out)
  th <- fit$estimates$theta
  out$sderr <- exp(th[1]); out$ar <- th[2] / sqrt(1 + th[2]^2)
  out$b1 <- fit$estimates$beta[if (D$fam == "cumulative") 1L else 2L]
  out$code <- fit$opt$convergence
  out$ll <- as.numeric(logLik(fit))
  info <- row_info(fit)
  out$info_med <- median(info); out$info_mean <- mean(info)
  out$info_q25 <- unname(quantile(info, 0.25))
  out$info_lo <- mean(info < 0.5)
  bk <- ns$frame_latres(fit$frame)[[1L]]
  env <- fit$obj$env
  H <- env$spHess(env$last.par.best, random = TRUE)
  # the fraction of a latent's prior variance the data leave: the Laplace
  # approximation's own account of how much of each latent the rows pin down
  out$n <- nrow(d)
  extra <- switch(D$fam,
    negbinomial = , gamma = fit$opt$par[names(fit$opt$par) == "betad"],
    NULL)
  ths <- if (D$fam == "cumulative") {
    as.numeric(ns$brms_ord_thresholds %||% NULL)
  }
  st <- switch(D$fam,
    cumulative = c(D$b1, log(D$s), D$ar / sqrt(1 - D$ar^2), cuts),
    negbinomial = , gamma = c(D$b1, log(D$s), D$ar / sqrt(1 - D$ar^2), D$b0, log(D$shape)),
    c(D$b1, log(D$s), D$ar / sqrt(1 - D$ar^2), D$b0))
  lo <- c(-Inf, log(0.02), -8, rep(-Inf, length(st) - 3L))
  hi <- c(Inf, log(10), 8, rep(Inf, length(st) - 3L))
  ex <- tryCatch(nlminb(st, exact_nll, d = d, lower = lo, upper = hi,
                        control = list(iter.max = 120, eval.max = 300)),
                 error = function(e) NULL)
  if (!is.null(ex)) {
    out$x_sderr <- exp(ex$par[2])
    out$x_ar <- ex$par[3] / sqrt(1 + ex$par[3]^2)
    out$x_b1 <- ex$par[1]
    out$x_ll <- -ex$objective
  }
  out
}
cols <- c("design", "seed", "sderr", "ar", "b1", "code", "ll", "info_med", "info_mean",
          "info_q25", "info_lo", "n", "x_sderr", "x_ar", "x_b1", "x_ll")
first <- !file.exists(outf)
for (sd in s1:s2) {
  o <- one(sd)
  for (cn in setdiff(cols, names(o))) o[[cn]] <- NA
  write.table(as.data.frame(o[cols]), outf, sep = "\t", row.names = FALSE,
              quote = FALSE, col.names = first, append = !first)
  first <- FALSE
}
