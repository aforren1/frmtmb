# Items 5620 and 5615: frmtmb_latres_laplace fires on lognormal, where
# the Laplace approximation is exact (its log density is quadratic in the
# linear predictor), and is decided by rounding between an unread mi()
# response and its weight-0 spelling. Data of dev/fitedge-mi-latres.R
# (seeds 40401 to 40405) and dev/review-fitedge-mi.R (lognormal ar(t, g),
# set.seed(73000 + seed)); plus the guard-present cases: poisson and
# bernoulli latents, where the approximation is biased.
#   Rscript dev/quadseam2-latres.R <libpath|base>
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
fit_w <- function(...) {
  w <- character()
  f <- tryCatch(withCallingHandlers(suppressMessages(frm(...)),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    }), error = function(e) e)
  list(fit = f, warn = w)
}
lat <- function(r) {
  if (inherits(r$fit, "error")) return("ERROR")
  sprintf("ll %.6f latres %d other %s", as.numeric(logLik(r$fit)),
          sum(r$warn == "frmtmb_latres_laplace"),
          paste(setdiff(r$warn, "frmtmb_latres_laplace"), collapse = ","))
}
for (seed in 40401:40405) {
  set.seed(seed)
  G <- 10L; K <- 6L
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  e <- as.vector(replicate(G, arima.sim(list(ar = 0.6), K, sd = 0.5)))
  d$y <- exp(0.6 + 0.4 * d$x + e + rnorm(nrow(d), 0, 0.3))
  miss <- sort(sample(nrow(d), 8))
  d$y[miss] <- NA
  r <- fit_w(bf(y | mi() ~ x + ar(t, g, cov = TRUE)), family = lognormal(),
             data = d)
  dw <- d; dw$w <- ifelse(is.na(d$y), 0, 1); dw$y[is.na(dw$y)] <- 1
  rw <- fit_w(bf(y | weights(w) ~ x + ar(t, g, cov = TRUE)),
              family = lognormal(), data = dw)
  cat(sprintf("LN %d mi: %s | w0: %s\n", seed, lat(r), lat(rw)))
}
# guard present: a poisson latent residual, counts of mean about 1, and a
# bernoulli one, where the warning is the point
for (seed in 1:3) {
  set.seed(500 + seed)
  G <- 20L; K <- 10L
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  e <- as.vector(replicate(G, arima.sim(list(ar = 0.5), K, sd = 0.8)))
  d$yp <- rpois(nrow(d), exp(0 + 0.3 * d$x + e))
  d$yb <- rbinom(nrow(d), 1, plogis(0.3 * d$x + e))
  rp <- fit_w(bf(yp ~ x + ar(t, g, cov = TRUE)), family = poisson(), data = d)
  rb <- fit_w(bf(yb ~ x + ar(t, g, cov = TRUE)), family = bernoulli(),
              data = d)
  cat(sprintf("POIS %d: %s | BERN %d: %s\n", seed, lat(rp), seed, lat(rb)))
}
cat("END\n")
