# M5, nlf(zi ~ ...), the Laplace value on brms's program at chosen points.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
Sys.setenv(FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/stan-cache")
suppressPackageStartupMessages(library(rstan))
hel <- new.env(parent = asNamespace("frmtmb"))
for (f in c("helper-brms.R", "helper-brms-methods.R", "helper-warnings.R")) sys.source(file.path("tests/testthat", f), envir = hel)
set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:30)); d <- d[order(d$g, d$t), ]
n <- nrow(d); d$x <- rnorm(n)
e <- as.vector(replicate(30, arima.sim(list(ar = 0.5), 5, sd = 0.5)))
d$count <- rpois(n, exp(0.8 + 0.4 * d$x + e))
fz <- count ~ x + ar(t, g, cov = TRUE)
fb <- brms::bf(fz) + brms::nlf(zi ~ exp(a) / (1 + exp(a))) + brms::lf(a ~ 1)
fam <- brms::brmsfamily("zero_inflated_poisson")
pz <- hel$brms_flat_prior(fb, d, fam)
code <- brms::stancode(fb, d, family = fam, prior = pz)
sdat <- brms::standata(fb, d, family = fam, prior = pz)
m <- hel$brms_stan_model(code)
sz <- rstan::sampling(m, data = sdat, chains = 0, refresh = 0)
f <- suppressWarnings(frm(frmtmb::bf(fz) + frmtmb::nlf(zi ~ exp(a) / (1 + exp(a))) + frmtmb::lf(a ~ 1), family = zero_inflated_poisson(), data = d))
env <- f$obj$env
print(names(f$opt$par))
print(sz@model_pars)
for (p in list(list(b0 = .5, b = .2, s = .6, ar = .4, a = -1), list(b0 = 1, b = -.1, s = .3, ar = -.5, a = -.6))) {
  nb <- grep("^b_.*a$", sz@model_pars, value = TRUE)[1]
  pl <- list(b = array(p$b, 1), Intercept = p$b0 + p$b * mean(d$x), ar = array(p$ar, 1), zerr = rep(0, n), sderr = p$s)
  pl[[nb]] <- array(p$a, 1)
  u <- rstan::unconstrain_pars(sz, pl)
  pl1 <- pl; pl1$zerr <- rep(1, n)
  iz <- which(rstan::unconstrain_pars(sz, pl1) != u)
  gr <- function(z) { u[iz] <- z; as.numeric(rstan::grad_log_prob(sz, u, adjust_transform = FALSE)[iz]) }
  lp <- function(z) { u[iz] <- z; as.numeric(rstan::log_prob(sz, u, adjust_transform = FALSE)) }
  hs <- function(z, h = 1e-5) { H <- matrix(0, n, n); for (i in 1:n) { e <- rep(0, n); e[i] <- h; H[i, ] <- (gr(z + e) - gr(z - e)) / (2 * h) }; (H + t(H)) / 2 }
  z <- rep(0, n)
  for (it in 1:100) {
    g <- gr(z)
    if (max(abs(g)) < 1e-9) break
    st <- solve(-hs(z), g)
    a <- 1
    while (a > 1e-6 && lp(z + a * st) < lp(z)) a <- a / 2
    z <- z + a * st
  }
  lb <- lp(z) + n / 2 * log(2 * pi) - .5 * as.numeric(determinant(-hs(z))$modulus)
  par <- f$opt$par
  nm <- names(par)
  par[nm == "beta"][1:2] <- c(p$b0, p$b)
  par[setdiff(seq_along(par), which(nm == "theta"))][3] <- p$a
  par[nm == "theta"] <- c(log(p$s), p$ar / sqrt(1 - p$ar^2))
  env$last.par.best[env$random] <- 0
  lf <- -as.numeric(f$obj$fn(par))
  cat(sprintf("brms %.8f frm %.8f diff %.2e\n", lb, lf, lf - lb))
}
