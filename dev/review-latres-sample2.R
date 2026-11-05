# Attack 3: the post-link forms, frm() against brms's posterior on data from the
# model (the DGP of review-latres-sim2.R). Default priors, flat shape on none.
# Usage: Rscript review-latres-sample2.R <design> <seed>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
args <- commandArgs(TRUE)
design <- args[[1]]
seed <- as.integer(args[[2]])
suppressPackageStartupMessages(library(rstan))
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
D <- switch(design,
  gamma_F = list(ff = Gamma(link = "log"), fb = brms::brmsfamily("Gamma", "log"), dist = "gamma", inv = exp, shape = 4, b0 = .5, b1 = .4, s = .5, ar = .5),
  probit_F = list(ff = bernoulli("probit"), fb = brms::brmsfamily("bernoulli", "probit"), dist = "bern", inv = pnorm, b0 = .2, b1 = .5, s = .5, ar = .2))
set.seed(seed)
G <- 25; K <- 6
d <- expand.grid(t = seq_len(K), g = seq_len(G))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$y <- NA_real_
for (g in seq_len(G)) {
  idx <- which(d$g == g)
  repeat {
    z <- rnorm(K)
    lag <- c(0, z[-K])
    mu <- D$inv(D$b0 + D$b1 * d$x[idx] + D$s * z) + D$ar * D$s * lag
    ok <- if (D$dist == "bern") all(mu > 0 & mu < 1) else all(mu > 0)
    if (ok) break
  }
  d$y[idx] <- if (D$dist == "gamma") rgamma(K, D$shape, D$shape / mu) else rbinom(K, 1, mu)
}
d$g <- factor(d$g)
m <- msg_of(frm(y ~ x + ar(t, g), family = D$ff, data = d))
f <- m$val
cat("frm: code", f$opt$convergence, " ll", as.numeric(logLik(f)), " nw", length(m$warn), "\n")
th <- f$estimates$theta
cat("frm: b", f$estimates$beta, " sderr", exp(th[1]), " ar", 2 / (1 + exp(-th[2])) - 1, " par", f$opt$par, "\n")
fb <- brms::bf(y ~ x + ar(t, g))
code <- brms::stancode(fb, d, family = D$fb)
sdat <- brms::standata(fb, d, family = D$fb)
mod <- hel$brms_stan_model(code)
sf <- rstan::sampling(mod, data = sdat, chains = 4, iter = 2000, warmup = 1000, seed = 11, refresh = 0,
                      init = function() list(sderr = .2, zerr = rep(0, n), ar = array(0, 1), b = array(0, 1), Intercept = D$b0))
pn <- c("b_Intercept", "b", "ar", "sderr", if (D$dist == "gamma") "shape")
print(round(rstan::summary(sf, pars = pn)$summary[, c("mean", "50%", "sd", "Rhat", "n_eff")], 3))
cat("divergences", sum(sapply(rstan::get_sampler_params(sf, inc_warmup = FALSE), function(x) sum(x[, "divergent__"]))), "\n")
