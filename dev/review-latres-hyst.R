# Attack 3: post-link AR. Is obj$fn a function of its argument? The same outer
# point from several feasible inner starts, on seeds of the Gamma log, Beta and
# probit forms. A second mode, or a start that fails, changes the value.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
run <- function(seed, fam, fam_lab, y) {
  set.seed(seed)
  d <- expand.grid(t = 1:6, g = factor(1:10))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- rnorm(n)
  e <- as.vector(replicate(10, arima.sim(list(ar = .5), 6, sd = .5)))
  d$pos <- rgamma(n, 4, 4 / exp(.5 + .4 * d$x + e))
  d$un <- pmin(pmax(rbeta(n, 8 * plogis(.3 * d$x + e), 8 * (1 - plogis(.3 * d$x + e))), 1e-3), 1 - 1e-3)
  d$bin <- rbinom(n, 1, pnorm(.3 * d$x + e)); d$cnt <- rpois(n, exp(.5 + .4 * d$x + e))
  m <- msg_of(frm(as.formula(paste(y, "~ x + ar(t, g)")), family = fam, data = d))
  if (!is.null(m$err)) return(data.frame(fam = fam_lab, seed = seed, code = NA, ll = NA, nan_starts = NA, spread = NA, n_values = NA, nw = NA))
  f <- m$val
  env <- f$obj$env
  par <- f$opt$par
  r <- env$random
  bb <- f$estimates$b
  vals <- numeric(0)
  nan <- 0L
  for (i in 1:25) {
    set.seed(seed * 100 + i)
    env$last.par.best[r] <- if (i == 1) 0 else if (i == 2) bb else bb + rnorm(length(r), 0, .25 * (i %% 5))
    v <- tryCatch(suppressWarnings(f$obj$fn(par)), error = function(e) NaN)
    if (is.finite(v)) vals <- c(vals, as.numeric(v)) else nan <- nan + 1L
  }
  data.frame(fam = fam_lab, seed = seed, code = f$opt$convergence, ll = as.numeric(logLik(f)),
             nan_starts = nan, spread = if (length(vals)) max(vals) - min(vals) else NA,
             lowest_nll_minus_fit = if (length(vals)) min(vals) + as.numeric(logLik(f)) else NA,
             n_values = length(unique(round(vals, 4))), nw = length(m$warn))
}
out <- list()
for (s in 1:12) {
  out[[length(out) + 1]] <- run(s, Gamma(link = "log"), "Gamma log", "pos")
  out[[length(out) + 1]] <- run(s, Beta(), "Beta logit", "un")
  out[[length(out) + 1]] <- run(s, bernoulli("probit"), "probit", "bin")
  out[[length(out) + 1]] <- run(s, poisson(), "poisson log control", "cnt")
}
o <- do.call(rbind, lapply(out, function(z) {
  for (cn in c("lowest_nll_minus_fit")) if (is.null(z[[cn]])) z[[cn]] <- NA
  z[c("fam", "seed", "code", "ll", "nan_starts", "spread", "lowest_nll_minus_fit", "n_values", "nw")]
}))
o <- o[order(o$fam, o$seed), ]
print(o, digits = 4, row.names = FALSE)
