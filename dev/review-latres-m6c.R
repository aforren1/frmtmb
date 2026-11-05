# Short check, item 4, brms's side, and M5's zi by the Laplace value. Same data as
# review-latres-m6b.R. brms's program with a held point, so epred reads the same
# parameters as frm()'s estimates, and with sampling.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages(library(rstan))
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
set.seed(167)
d <- expand.grid(t = 1:6, g = factor(1:25))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(25, arima.sim(list(ar = .5), 6, sd = .5)))
d$cnt <- rpois(n, exp(.8 + .3 * d$x + e))
f <- frm(cnt ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d)
b <- f$estimates$beta
th <- f$estimates$theta
nd3 <- data.frame(t = 1:6, g = factor(1), x = 0)
pri <- c(brms::set_prior("", class = "sderr"), brms::set_prior("", class = "Intercept"))
fb <- brms::bf(cnt ~ x + ar(t, g, cov = TRUE))
fam <- brms::brmsfamily("poisson")
code <- brms::stancode(fb, d, family = fam, prior = pri)
sdat <- brms::standata(fb, d, family = fam, prior = pri)
mod <- hel$brms_stan_model(code)
pp <- list(b = array(b[2], 1), Intercept = b[1] + b[2] * mean(d$x), ar = array(th[2] / sqrt(1 + th[2]^2), 1), zerr = rep(0, n), sderr = exp(th[1]))
sf <- rstan::sampling(mod, data = sdat, algorithm = "Fixed_param", chains = 1, iter = 8000, warmup = 0, init = list(pp), refresh = 0, seed = 1)
scaf <- brms::brm(fb, d, family = fam, prior = pri, empty = TRUE)
scaf$fit <- sf
bh <- brms::rename_pars(scaf)
set.seed(9)
ep <- brms::posterior_epred(bh, newdata = nd3)
cat("held point: brms epred mean", mean(ep), " MC se", sd(rowMeans(ep)) / sqrt(nrow(ep)), " frm fitted", mean(fitted(f, newdata = nd3)[, 1]), " marginaleffects", mean(marginaleffects::predictions(f, newdata = nd3, vcov = FALSE)$estimate), "\n")
mh <- msg_of(marginaleffects::predictions(bh, newdata = nd3))
cat("marginaleffects on the held brms fit:", if (is.null(mh$err)) round(mean(mh$val$estimate), 4) else substr(mh$err, 1, 200), "\n")
cat("  (a median over draws, with err drawn per draw:", round(median(ep[, 1]), 4), " mean:", round(mean(ep[, 1]), 4), ")\n")
s1 <- msg_of(marginaleffects::avg_slopes(bh, variables = "x"))
s2 <- marginaleffects::avg_slopes(f, variables = "x", vcov = FALSE)
cat("avg_slopes at the held point: brms", if (is.null(s1$err)) s1$val$estimate else substr(s1$err, 1, 160), " frm", s2$estimate, "\n")
# the slope brms's mean over draws would read: b1 * mean(E)
ep2 <- brms::posterior_epred(bh, newdata = d)
cat("mean of brms epred on d", mean(ep2), " frm", mean(fitted(f, newdata = d)[, 1]), "; b1 x that:", b[2] * mean(ep2), "\n")

cat("\n== M5, nlf(zi)\n")
set.seed(167)
d2 <- expand.grid(t = 1:4, g = factor(1:8))
d2 <- d2[order(d2$g, d2$t), ]
m <- nrow(d2)
d2$x <- rnorm(m)
d2$z <- rnbinom(m, size = 3, mu = exp(1 + .3 * d2$x)) * rbinom(m, 1, .75)
pz <- c(brms::set_prior("", class = "sderr"), brms::set_prior("", class = "Intercept"))
fz <- brms::bf(z ~ x + ar(t, g, cov = TRUE), brms::nlf(zi ~ exp(a) / (1 + exp(a))), a ~ 1)
famz <- brms::brmsfamily("zero_inflated_poisson")
code <- brms::stancode(fz, d2, family = famz, prior = pz)
sdat <- brms::standata(fz, d2, family = famz, prior = pz)
mz <- hel$brms_stan_model(code)
sz <- rstan::sampling(mz, data = sdat, chains = 0, refresh = 0)
ffz <- suppressWarnings(frm(frmtmb::bf(z ~ x + ar(t, g, cov = TRUE)) + frmtmb::nlf(zi ~ exp(a) / (1 + exp(a))) + frmtmb::lf(a ~ 1), family = zero_inflated_poisson(), data = d2))
env <- ffz$obj$env
for (p in list(list(b0 = .5, b = .2, s = .6, ar = .4, a = -1), list(b0 = 1, b = -.1, s = .3, ar = -.5, a = -.6))) {
  pl <- list(b = array(p$b, 1), Intercept = p$b0 + p$b * mean(d2$x), ar = array(p$ar, 1), zerr = rep(0, m), sderr = p$s, b_zi_a = array(p$a, 1))
  if (!grepl("b_zi_a", paste(names(sz@par_dims), collapse = " "))) {
    pl$b_a <- pl$b_zi_a
    pl$b_zi_a <- NULL
  }
  u <- rstan::unconstrain_pars(sz, pl)
  pl1 <- pl
  pl1$zerr <- rep(1, m)
  iz <- which(rstan::unconstrain_pars(sz, pl1) != u)
  gr <- function(z) { u[iz] <- z; as.numeric(rstan::grad_log_prob(sz, u, adjust_transform = FALSE)[iz]) }
  lp <- function(z) { u[iz] <- z; as.numeric(rstan::log_prob(sz, u, adjust_transform = FALSE)) }
  hs <- function(z, h = 1e-5) { H <- matrix(0, m, m); for (i in 1:m) { e <- rep(0, m); e[i] <- h; H[i, ] <- (gr(z + e) - gr(z - e)) / (2 * h) }; (H + t(H)) / 2 }
  z <- rep(0, m)
  for (it in 1:100) {
    g <- gr(z)
    if (max(abs(g)) < 1e-9) break
    st <- solve(-hs(z), g)
    a <- 1
    while (a > 1e-6 && lp(z + a * st) < lp(z)) a <- a / 2
    z <- z + a * st
  }
  lb <- lp(z) + m / 2 * log(2 * pi) - .5 * as.numeric(determinant(-hs(z))$modulus)
  par <- ffz$opt$par
  nm <- names(par)
  par[nm == "beta"][1:2] <- c(p$b0, p$b)
  par[setdiff(seq_along(par), which(nm == "theta"))][3] <- p$a
  par[nm == "theta"] <- c(log(p$s), p$ar / sqrt(1 - p$ar^2))
  env$last.par.best[env$random] <- 0
  lf <- -as.numeric(ffz$obj$fn(par))
  cat(sprintf("brms %.8f frm %.8f diff %.2e\n", lb, lf, lf - lb))
}
