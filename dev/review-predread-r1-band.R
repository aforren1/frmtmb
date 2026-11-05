# Short check: fitted()'s interval on a me() row when more than one predictor reads the
# term. Run on the lane's build and on mutant band_var, the proposed one-line change.
# Seed 8401, the data of review-predread-r1-check7.R.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(8401)
n <- 500
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4)
d$xo <- d$lat + rnorm(n, 0, 0.4)
eta <- 0.4 + 0.7 * d$lat
d$cnt <- rpois(n, exp(eta))
d$zc <- ifelse(runif(n) < plogis(-1 + 0.8 * d$lat), 0L, rpois(n, 3))
comp <- rbinom(n, 1, 0.5)
d$mixc <- rpois(n, ifelse(comp == 1, exp(1.5 + 0.4 * d$lat), exp(0.2 + 0.3 * d$lat)))
nd <- data.frame(xo = c(1, 1), sx = c(0.4, 1.2), z = 0.3)
fmt <- function(v) paste(sprintf("%.3f", v), collapse = " ")
for (sp in list(list("bf(cnt ~ me(xo, sx))", poisson()),
                list("bf(zc ~ me(xo, sx), zi ~ me(xo, sx))", zero_inflated_poisson()),
                list("bf(mixc ~ me(xo, sx))", mixture(poisson(), poisson())))) {
  fit <- quiet(frm(rev_ev(sp[[1]], "frmtmb"), family = sp[[2]], data = d))
  fv <- quiet(fitted(fit, newdata = nd))
  ds <- rev_draws_at(fit, 8000L)
  set.seed(9)
  fe <- quiet(frmtmb.sample::posterior_epred(ds, newdata = nd))
  cat(sp[[1]], "\n  fitted Q ", fmt(fv[, "Q2.5"]), "/", fmt(fv[, "Q97.5"]), " est ", fmt(fv[, "Estimate"]), " se ", fmt(fv[, "Est.Error"]),
      "\n  draws  q ", fmt(apply(fe, 2, quantile, 0.025)), "/", fmt(apply(fe, 2, quantile, 0.975)), " sd ", fmt(apply(fe, 2, sd)), "\n")
}
