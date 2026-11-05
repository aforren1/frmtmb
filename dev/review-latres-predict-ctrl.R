# predict()'s blow-up: the same model written as the latent term and as an
# ar1() random effect, same data, same seeds, and predict() with and without the
# parameter draw. The control tells whether the lane's route is the cause.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
res <- list()
for (seed in c(1, 3, 7, 9, 13, 17, 19)) {
  set.seed(seed)
  d <- expand.grid(t = 1:6, g = factor(1:10))
  n <- nrow(d)
  d$x <- rnorm(n)
  d$ft <- factor(d$t)
  e <- as.vector(replicate(10, arima.sim(list(ar = .4), 6, sd = .25)))
  d$y <- rpois(n, exp(0.5 + 0.3 * d$x + e))
  fl <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
  fr <- suppressWarnings(frm(y ~ x + ar1(ft + 0 | g), family = poisson(), data = d))
  pl <- suppressWarnings(predict(fl, ndraws = 500))[, 1]
  pr <- suppressWarnings(predict(fr, ndraws = 500))[, 1]
  pl0 <- suppressWarnings(predict(fl, ndraws = 500, propagate_error = FALSE))[, 1]
  vl <- sqrt(diag(as.matrix(vcov(fl, full = TRUE))))
  vr <- sqrt(diag(as.matrix(vcov(fr, full = TRUE))))
  res[[length(res) + 1]] <- data.frame(seed = seed, ll_lat = as.numeric(logLik(fl)), ll_ar1 = as.numeric(logLik(fr)),
    sd_lat = exp(fl$estimates$theta[1]), se_logsd_lat = vl["theta_1"], se_logsd_ar1 = vr["theta_1"],
    max_y = max(d$y), max_pred_lat = max(pl), max_pred_ar1 = max(pr), max_pred_lat_noprop = max(pl0))
}
print(do.call(rbind, res), digits = 4, row.names = FALSE)
cat("\n`predict(fit, newdata)` on the same seeds, a row of mean exp(eta + err):\n")
