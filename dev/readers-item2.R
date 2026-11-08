# Item 2: conditional_effects()'s band for a category under me(), against
# fitted()'s Q columns at the same row and against the draw's own quantiles at
# the estimates (brms's posterior_epred at fixed values, by Monte Carlo, 2e4
# draws, seed 11). Data: the design of dev/review-predread-r1-check1.R, seed
# 7971; row xo = 1, sx = 1.2.
source("dev/readers-setup.R")
q <- function(e) tryCatch(suppressWarnings(e), error = function(err) {
  cat("ERROR:", conditionMessage(err), "\n"); NULL
})
set.seed(7971)
n <- 400
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4)
d$xo <- d$lat + rnorm(n, 0, 0.4)
eta <- 0.4 + 0.7 * d$lat
d$ord <- as.integer(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
                        labels = FALSE))
d$cat <- as.integer(cut(eta + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
nd <- data.frame(xo = 1, sx = 1.2)
for (sp in list(list(bf(ord ~ me(xo, sx)), cumulative()),
                list(bf(cat ~ me(xo, sx)), categorical()))) {
  f <- frm(sp[[1]], family = sp[[2]], data = d)
  fv <- fitted(f, newdata = nd)
  ce <- q(conditional_effects(f, effects = "xo", categorical = TRUE,
                              int_conditions = list(xo = 1),
                              conditions = data.frame(sx = 1.2)))
  # the draw at the estimates: the predictor moves by b * (X - xo)
  set.seed(11)
  zz <- rnorm(2e5)
  K <- dim(fv)[3L]
  cat("\n===", sp[[2]]$family, "\n")
  # draws of the latent value, read through the fit at X = xo + sx * z
  Xs <- nd$xo + nd$sx * zz[1:20000]
  Pm <- as.matrix(frm_linpred(f, newdata = data.frame(xo = Xs, sx = 1.2),
                              type = "response"))
  for (k in seq_len(K)) {
    cat(sprintf(paste("cat %d  fitted Q %.3f to %.3f | ce %s | draw at the",
                      "estimate %.3f to %.3f\n"), k,
                fv[1, "Q2.5", k], fv[1, "Q97.5", k],
                if (is.null(ce)) "ERROR" else
                  sprintf("%.3f to %.3f", ce[[1]]$lower__[k],
                          ce[[1]]$upper__[k]),
                quantile(Pm[, k], 0.025), quantile(Pm[, k], 0.975)))
  }
}
