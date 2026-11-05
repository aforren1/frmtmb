# Short check: M1 and M3 on routes the family table does not reach: me() in a body, in a
# dpar the mean reads, beside trunc(), in a mixture, beside cs(). fitted() against frmtmb's
# own draws (the D1 reading) and brms's at fixed values, x = 1, sdx 1.2. Seed 8401.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(8401)
n <- 500
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4)
d$xo <- d$lat + rnorm(n, 0, 0.4)
eta <- 0.4 + 0.7 * d$lat
d$cnt <- rpois(n, exp(eta))
d$pos <- exp(0.3 + 0.2 * d$z + rnorm(n, 0, exp(-0.8 + 0.5 * d$lat)))
d$zc <- ifelse(runif(n) < plogis(-1 + 0.8 * d$lat), 0L, rpois(n, 3))
d$tr <- pmin(d$cnt, 12L)
comp <- rbinom(n, 1, 0.5)
d$mixc <- rpois(n, ifelse(comp == 1, exp(1.5 + 0.4 * d$lat), exp(0.2 + 0.3 * d$lat)))
d$ord <- as.integer(cut(0.7 * d$lat + ifelse(d$z > 0, 1.2, -1.2) * rlogis(n), c(-Inf, -0.5, 0.8, 2, Inf), labels = FALSE))
nd <- data.frame(xo = c(1, 1), sx = c(0.4, 1.2), z = 0.3)
specs <- list(
  body_poisson = c("bf(cnt ~ a + b * z, a ~ me(xo, sx), b ~ 1, nl = TRUE)", "poisson()"),
  body_exp = c("bf(cnt ~ a + exp(b) * z, a ~ me(xo, sx), b ~ 1, nl = TRUE)", "poisson()"),
  lognormal_sigma = c("bf(pos ~ z, sigma ~ me(xo, sx))", "lognormal()"),
  zi_me = c("bf(zc ~ z, zi ~ me(xo, sx))", "zero_inflated_poisson()"),
  mu_and_zi = c("bf(zc ~ me(xo, sx), zi ~ me(xo, sx))", "zero_inflated_poisson()"),
  trunc = c("bf(tr | trunc(ub = 12) ~ me(xo, sx))", "poisson()"),
  mixture = c("bf(mixc ~ me(xo, sx))", "mixture(poisson(), poisson())"),
  cs_me = c("bf(ord ~ me(xo, sx) + cs(z))", "acat()"))
fmt <- function(v) paste(sprintf("%.4f", v), collapse = " ")
one <- function(a) if (length(dim(a)) == 3L) a[, , 1] else a
for (nm in names(specs)) {
  cat("\n===", nm, specs[[nm]], "\n")
  fit <- tryCatch(quiet(frm(rev_ev(specs[[nm]][1], "frmtmb"), family = rev_ev(specs[[nm]][2], "frmtmb"), data = d)), error = function(e) e)
  if (inherits(fit, "error")) {
    cat("frm ERR", substr(conditionMessage(fit), 1, 170), "\n")
    next
  }
  fv <- tryCatch(quiet(fitted(fit, newdata = nd)), error = function(e) e)
  if (inherits(fv, "error")) {
    cat("fitted ERR", substr(conditionMessage(fv), 1, 170), "\n")
  } else {
    if (length(dim(fv)) == 3L) fv <- fv[, , 1]
    cat("fitted      ", fmt(fv[, "Estimate"]), " se ", fmt(fv[, "Est.Error"]), " Q ", fmt(fv[, "Q2.5"]), "/", fmt(fv[, "Q97.5"]), "\n")
  }
  ds <- rev_draws_at(fit, 8000L)
  set.seed(9)
  fe <- tryCatch(one(quiet(frmtmb.sample::posterior_epred(ds, newdata = nd))), error = function(e) e)
  if (inherits(fe, "error")) {
    cat("draws ERR", substr(conditionMessage(fe), 1, 170), "\n")
  } else {
    cat("draws epred ", fmt(colMeans(fe)), " sd ", fmt(apply(fe, 2, sd)), " q ", fmt(apply(fe, 2, quantile, 0.025)), "/", fmt(apply(fe, 2, quantile, 0.975)), "\n")
  }
  fx <- tryCatch(quiet(rev_E$brms_fixed_fit(rev_ev(specs[[nm]][1], "brms"), rev_ev(specs[[nm]][2], "brms"), d, fit, ndraws = 8000)), error = function(e) e)
  if (inherits(fx, "error")) {
    cat("brms ERR", substr(conditionMessage(fx), 1, 120), "\n")
  } else {
    be <- one(quiet(posterior_epred(fx, newdata = nd)))
    cat("brms epred  ", fmt(colMeans(be)), " sd ", fmt(apply(be, 2, sd)), " q ", fmt(apply(be, 2, quantile, 0.025)), "/", fmt(apply(be, 2, quantile, 0.975)), "\n")
  }
}
