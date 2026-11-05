# Attack 2, the boundary. (a) a random-intercept DGP, where ar -> 1 and
# sderr -> 0 with the marginal sd finite. (b) no latent at all, sderr -> 0.
# What does the user see from frm(), summary(), predict(), simulate()?
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
run <- function(seed, dgp, term = "ar(t, g, cov = TRUE)", fam = poisson()) {
  set.seed(seed)
  d <- expand.grid(t = 1:6, g = factor(1:10))
  n <- nrow(d)
  d$x <- rnorm(n)
  re <- switch(dgp, ri = rep(rnorm(10, 0, .5), each = 6), none = 0)
  eta <- 0.5 + 0.3 * d$x + re
  d$y <- if (fam$family == "poisson") rpois(n, exp(eta)) else rbinom(n, 1, plogis(eta))
  m <- msg_of(frm(as.formula(paste("y ~ x +", term)), family = fam, data = d))
  f <- m$val
  if (is.null(f)) return(data.frame(seed = seed, dgp = dgp, err = substr(m$err, 1, 100)))
  th <- f$estimates$theta
  sm <- summary(f)$cor_pars
  ms <- msg_of(capture.output(p <- predict(f, ndraws = 500)))
  pm <- msg_of(predict(f, ndraws = 500))
  sim <- msg_of(simulate(f, nsim = 200))
  mx <- if (!is.null(pm$val)) max(pm$val[, 1]) else NA
  simmax <- if (!is.null(sim$val)) max(as.matrix(sim$val), na.rm = TRUE) else NA
  sd_th <- tryCatch(sqrt(diag(as.matrix(vcov(f, full = TRUE))))[c("theta_1", "theta_2")], error = function(e) c(NA, NA))
  bm <- tryCatch(isTRUE(f$cache$boundary), error = function(e) NA)
  data.frame(seed = seed, dgp = dgp, code = f$opt$convergence, ll = as.numeric(logLik(f)),
             th1 = th[1], th2 = th[2], sderr = exp(th[1]), ar = sm[1, 1],
             se_ar = sm[1, 2], lo_ar = sm[1, 3], hi_ar = sm[1, 4], se_ls = sd_th[1], se_t = sd_th[2],
             max_y = max(d$y), pred_max = mx, sim_max = simmax,
             nw = length(m$warn), npw = length(pm$warn),
             w1 = substr(m$warn[1], 1, 70), diag_clean = tryCatch(length(capture.output(print(diagnose(f)))), error = function(e) NA))
}
ri <- do.call(rbind, lapply(1:15, run, dgp = "ri"))
print(ri[, c("seed", "code", "ll", "sderr", "ar", "se_ar", "lo_ar", "hi_ar", "se_ls", "se_t", "max_y", "pred_max", "sim_max", "nw", "npw")], digits = 4, row.names = FALSE)
print(table(substr(ri$w1, 1, 70), useNA = "ifany"))
none <- do.call(rbind, lapply(1:15, run, dgp = "none"))
print(none[, c("seed", "code", "ll", "sderr", "ar", "se_ar", "se_ls", "se_t", "max_y", "pred_max", "sim_max", "nw", "npw")], digits = 4, row.names = FALSE)
print(table(substr(none$w1, 1, 70), useNA = "ifany"))
# the first fit of the review, as a user reads it
set.seed(1)
d <- expand.grid(t = 1:6, g = factor(1:10))
d$x <- rnorm(nrow(d))
d$y <- rpois(nrow(d), exp(0.5 + 0.3 * d$x + rep(rnorm(10, 0, .5), each = 6)))
m <- msg_of(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
cat("\nwarnings:\n"); cat(m$warn, sep = "\n")
print(summary(m$val))
print(diagnose(m$val)$singular)
pp <- msg_of(predict(m$val, ndraws = 1000))
print(head(pp$val, 3)); cat(substr(pp$warn, 1, 300), sep = "\n")
pp2 <- msg_of(predict(m$val, newdata = d[1:3, ], ndraws = 1000))
print(pp2$val); cat(substr(pp2$warn, 1, 300), sep = "\n")
f2 <- frm(y ~ x + cosy(t, g), family = poisson(), data = d)
cat("cosy logLik", as.numeric(logLik(f2)), " ar logLik", as.numeric(logLik(m$val)), " (1|g) logLik", as.numeric(logLik(frm(y ~ x + (1 | g), family = poisson(), data = d))), "\n")
print(summary(f2)$cor_pars)
