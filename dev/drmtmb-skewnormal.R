# Root cause of the skew-normal logLik gap between frmtmb and drmTMB.
# frmtmb's alpha start takes the sign of the RAW response's sample
# skewness. When covariates or group effects make the marginal skew
# opposite to the residual skew, the start lands on the wrong side of
# alpha = 0, which is a stationary point in the mean parameterization
# (dgamma1/dalpha = 0 there), and nlminb stalls at it.
source("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-agree-lib.R")
skew <- function(v) mean((v - mean(v))^3) / sd(v)^3
start_alpha <- function(y) { m3 <- skew(y); 2 * sign(m3) + 0.5 * m3 }

d <- sim_grouped(101)
e <- d$ysn - 1 - 0.5 * d$x
cat("seed 101: skew(raw y) =", skew(d$ysn), " skew(y - true mean) =",
    skew(e), " frmtmb alpha start =", start_alpha(d$ysn), "\n")
ff <- frm(bf(ysn ~ x, sigma ~ 1, alpha ~ 1), family = skew_normal(), data = d)
fd <- drmTMB::drmTMB(dbf(ysn ~ x, sigma ~ 1, nu ~ 1),
                     family = drmTMB::skew_normal(), data = d)
ff_s <- frm(bf(ysn ~ x, sigma ~ 1, alpha ~ 1), family = skew_normal(),
            data = d, start = list(betad = c(-0.3, 1)))
op <- options(digits = 12)
print(rbind(frm_default = c(ll = logLik(ff), ff$opt$par),
            frm_start_pos = c(ll = logLik(ff_s), ff_s$opt$par),
            drm = c(ll = logLik(fd), fd$opt$par)))
options(op)

# Profile of the shared objective along alpha, other parameters at
# their optimum for each alpha, from frmtmb's own objective.
prof <- sapply(c(-1, -0.5, -0.1, 0, 0.1, 0.3, 0.36, 0.5, 1), function(a) {
  f <- frm(bf(ysn ~ x, sigma ~ 1, alpha ~ 1), family = skew_normal(),
           data = d, start = list(betad = c(-0.3, a)))
  # fix alpha by evaluating the profile with an inner optimize over the
  # other three parameters
  fn <- function(p) f$obj$fn(c(p[1:3], a))
  o <- nlminb(ff$opt$par[1:3], fn)
  c(alpha = a, loglik = -o$objective)
})
print(t(prof), digits = 10)

# Rate over seeds: a minimal design, y = x + skew-normal residual with
# positive skew, and a covariate large enough to flip the marginal skew.
rate <- t(sapply(1:40, function(s) {
  set.seed(s)
  n <- 200
  # Negative skew in the covariate, positive in the residual. An earlier
  # version drew an unused rnorm(n) here, which shifted the stream and
  # gave 33 of 40 rather than the 28 of 40 this design gives;
  # dev/drmtmb-skewstart.R measures both streams.
  xs <- -abs(rnorm(n)) * 3
  y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  dd <- data.frame(y = y, xs = xs)
  f1 <- frm(bf(y ~ xs, sigma ~ 1, alpha ~ 1), family = skew_normal(),
            data = dd)
  f2 <- drmTMB::drmTMB(dbf(y ~ xs, sigma ~ 1, nu ~ 1),
                       family = drmTMB::skew_normal(), data = dd)
  c(seed = s, skew_raw = skew(y), skew_resid = skew(resid(lm(y ~ xs))),
    ll_gap = as.numeric(logLik(f1)) - as.numeric(logLik(f2)),
    alpha_frm = f1$opt$par[[4]], nu_drm = f2$opt$par[[4]])
}))
print(round(rate, 5))
cat("seeds with frmtmb below drmTMB by more than 1e-6:",
    sum(rate[, "ll_gap"] < -1e-6), "of", nrow(rate), "\n")
cat("of those, raw skew sign opposite residual skew sign:",
    sum(rate[, "ll_gap"] < -1e-6 &
          sign(rate[, "skew_raw"]) != sign(rate[, "skew_resid"])), "\n")
cat("seeds with opposite signs:",
    sum(sign(rate[, "skew_raw"]) != sign(rate[, "skew_resid"])), "\n")
