# Attack 2: Laplace against the exact marginal likelihood, and SE coverage.
#
# The exact reference is the forward algorithm on a fine grid: the latent AR(1)
# is one-dimensional and Markov, so for each group the marginal likelihood is a
# chain of one-dimensional integrals, done by the trapezoid rule on a grid of
# `ngrid` points over +-6 marginal sds (spectral accuracy for smooth integrands).
# This takes no Laplace step anywhere, and reads only this file's model.
#
# Usage: Rscript review-latres-sim.R <design> <seed_from> <seed_to> <outfile>
# Each dataset: frm() (estimates, SE, logLik, code, warnings) and the exact MLE
# by nlminb on the forward likelihood.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
args <- commandArgs(TRUE)
design <- args[[1]]
s1 <- as.integer(args[[2]])
s2 <- as.integer(args[[3]])
outf <- args[[4]]

designs <- list(
  # small counts: mean about 1
  pois = list(fam = poisson(), G = 10, K = 15, b0 = -0.2, b1 = 0.5, s = 0.5, ar = 0.6),
  # strong AR, small counts
  pois_ar9 = list(fam = poisson(), G = 10, K = 15, b0 = -0.2, b1 = 0.5, s = 0.3, ar = 0.9),
  # bigger sd on small counts
  pois_s1 = list(fam = poisson(), G = 10, K = 15, b0 = -0.5, b1 = 0.5, s = 1, ar = 0.6),
  # binary, moderate
  bern = list(fam = bernoulli(), G = 20, K = 10, b0 = 0.2, b1 = 0.5, s = 1, ar = 0.6),
  # binary long series
  bern_long = list(fam = bernoulli(), G = 5, K = 60, b0 = 0.2, b1 = 0.5, s = 1, ar = 0.8),
  # binomial, few trials
  binom3 = list(fam = binomial(), G = 15, K = 10, b0 = 0.2, b1 = 0.5, s = 0.8, ar = 0.6, tr = 3),
  # large counts: where Laplace should be good
  pois_big = list(fam = poisson(), G = 10, K = 15, b0 = 3, b1 = 0.3, s = 0.5, ar = 0.6))
D <- designs[[design]]

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
  if (D$fam$family == "poisson") {
    d$y <- rpois(nrow(d), exp(eta))
  } else if (D$fam$family == "bernoulli") {
    d$y <- rbinom(nrow(d), 1, plogis(eta))
  } else {
    d$y <- rbinom(nrow(d), D$tr, plogis(eta))
    d$tr <- D$tr
  }
  d$g <- factor(d$g)
  d
}

# log p(y | eta) on the grid: one row per obs, returned per grid point
row_dens <- function(y, eta) {
  if (D$fam$family == "poisson") return(dpois(y, exp(eta)))
  if (D$fam$family == "bernoulli") return(dbinom(y, 1, plogis(eta)))
  dbinom(y, D$tr, plogis(eta))
}

exact_nll <- function(par, d, ngrid = 121) {
  b0 <- par[1]; b1 <- par[2]; s <- exp(par[3])
  th <- par[4]
  ar <- th / sqrt(1 + th^2)
  msd <- s * sqrt(1 + th^2)      # s / sqrt(1 - ar^2)
  u <- seq(-6.5, 6.5, length.out = ngrid)
  h <- u[2] - u[1]
  grid <- u * msd
  # transition density on the grid, with the proper mass for a sharp kernel:
  # when s is small against the spacing the rule loses accuracy, so the grid
  # is refined by the caller's ngrid
  Tm <- outer(grid, grid, function(a, b) dnorm(b, ar * a, s)) * (h * msd)
  ll <- 0
  for (idx in split(seq_len(nrow(d)), d$g)) {
    f <- dnorm(u) * h
    for (k in seq_along(idx)) {
      i <- idx[k]
      if (k > 1L) f <- as.vector(crossprod(f, Tm))
      f <- f * row_dens(d$y[i], b0 + b1 * d$x[i] + grid)
      c0 <- sum(f)
      if (!is.finite(c0) || c0 <= 0) return(1e10)
      ll <- ll + log(c0)
      f <- f / c0
    }
  }
  -ll
}

fit_one <- function(seed) {
  d <- simd(seed)
  lhs <- if (D$fam$family == "binomial") "y | trials(tr)" else "y"
  form <- as.formula(paste(lhs, "~ x + ar(t, g, cov = TRUE)"))
  t0 <- proc.time()[[3]]
  m <- msg_of(frm(form, data = d, family = D$fam))
  tfit <- proc.time()[[3]] - t0
  out <- data.frame(design = design, seed = seed, tfit = tfit)
  if (!is.null(m$err)) {
    out$err <- substr(m$err, 1, 200)
    return(out)
  }
  f <- m$val
  th <- f$estimates$theta
  out$err <- ""
  out$b0 <- f$estimates$beta[1]
  out$b1 <- f$estimates$beta[2]
  out$ls <- th[1]
  out$t1 <- th[2]
  out$ar <- th[2] / sqrt(1 + th[2]^2)
  out$sderr <- exp(th[1])
  out$ll <- as.numeric(logLik(f))
  out$code <- f$opt$convergence
  vc <- tryCatch(as.matrix(vcov(f, full = TRUE)), error = function(e) NULL)
  se <- tryCatch(sqrt(diag(vc)), error = function(e) rep(NA, 4))
  out$se_b0 <- se[1]
  out$se_b1 <- se[2]
  out$se_ls <- se[3]
  out$se_t1 <- se[4]
  rs <- tryCatch(summary(f)$cor_pars, error = function(e) NULL)
  ci <- tryCatch({
    suppressWarnings(confint_varcorr(f))
  }, error = function(e) NULL)
  if (is.data.frame(ci)) {
    getci <- function(col, tm) {
      i <- which(ci$term == tm)
      if (length(i)) ci[[col]][i[1]] else NA_real_
    }
    out$sderr_lo <- getci("lwr", "sderr")
    out$sderr_hi <- getci("upr", "sderr")
    out$ar_lo <- getci("lwr", "ar[1]")
    out$ar_hi <- getci("upr", "ar[1]")
  }
  out$nwarn <- length(m$warn)
  out$warn <- substr(paste(m$warn, collapse = " | "), 1, 300)
  # the exact MLE from the truth and from frm()'s point
  st <- c(D$b0, D$b1, log(D$s), D$ar / sqrt(1 - D$ar^2))
  t0 <- proc.time()[[3]]
  ex <- tryCatch(nlminb(st, exact_nll, d = d,
                        lower = c(-Inf, -Inf, log(0.02), -8),
                        upper = c(Inf, Inf, log(10), 8),
                        control = list(iter.max = 200, eval.max = 400)),
                 error = function(e) NULL)
  if (!is.null(ex)) {
    out$tex <- proc.time()[[3]] - t0
    out$x_b0 <- ex$par[1]
    out$x_b1 <- ex$par[2]
    out$x_ls <- ex$par[3]
    out$x_t1 <- ex$par[4]
    out$x_ar <- ex$par[4] / sqrt(1 + ex$par[4]^2)
    out$x_sderr <- exp(ex$par[3])
    out$x_ll <- -ex$objective
    out$x_code <- ex$convergence
    # the exact value at frm()'s own point
    pf <- c(out$b0, out$b1, out$ls, pmin(pmax(out$t1, -8), 8))
    out$x_at_frm <- -exact_nll(pf, d)
  }
  out
}

first <- TRUE
for (sd in s1:s2) {
  o <- fit_one(sd)
  all_cols <- c("design","seed","tfit","err","b0","b1","ls","t1","ar","sderr","ll","code",
                "se_b0","se_b1","se_ls","se_t1","sderr_lo","sderr_hi","ar_lo","ar_hi",
                "nwarn","warn","tex","x_b0","x_b1","x_ls","x_t1","x_ar","x_sderr",
                "x_ll","x_code","x_at_frm")
  for (cn in setdiff(all_cols, names(o))) o[[cn]] <- NA
  o <- o[all_cols]
  write.table(o, outf, sep = "\t", row.names = FALSE, quote = TRUE,
              col.names = first, append = !first)
  first <- FALSE
}
