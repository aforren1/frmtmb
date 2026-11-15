# Review of lane quadseam2, attack 5: latres_quadratic(), the gate that
# drops frmtmb_latres_laplace where every row's log density is quadratic
# in eta. Is it scale-free, and does it drop the warning on a family that
# is not quadratic? Per fit: count of the warning, the gate's answer
# (lane only), and logLik.
#   Rscript quadseam2-review-latres.R <libpath|base>
a <- commandArgs(TRUE)
if (a[1] == "lib") a[1] <- "C:/Users/adf44/source/r/wt-quadseam2-lib"
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
has_q <- exists("latres_quadratic", ns)
fit_w <- function(...) {
  w <- character()
  f <- tryCatch(withCallingHandlers(suppressMessages(frm(...)),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    }), error = function(e) e)
  list(fit = f, warn = w)
}
lat <- function(r) {
  if (inherits(r$fit, "error")) {
    return(paste("ERROR", substr(conditionMessage(r$fit), 1, 90)))
  }
  q <- if (has_q) {
    tryCatch(ns$latres_quadratic(r$fit, ns$frame_latres(r$fit$frame)[[1L]]),
             error = function(e) NA)
  } else "-"
  sh <- tryCatch(ns$latres_shift_tk(r$fit, ns$frame_latres(r$fit$frame)[[1L]]),
                 error = function(e) NA)
  sprintf("ll %.6f latres %d quad %s shift %.3g other %s",
          as.numeric(logLik(r$fit)),
          sum(r$warn == "frmtmb_latres_laplace"), q, sh,
          paste(setdiff(r$warn, "frmtmb_latres_laplace"), collapse = ","))
}
ln_data <- function(seed) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(6L), g = factor(seq_len(10L)))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  e <- as.vector(replicate(10L, arima.sim(list(ar = 0.6), 6L, sd = 0.5)))
  d$y <- exp(0.6 + 0.4 * d$x + e + rnorm(nrow(d), 0, 0.3))
  d$y[sort(sample(nrow(d), 8))] <- NA
  d
}
for (sc in c(1e-8, 1e8)) for (seed in 40401:40403) {
  d <- ln_data(seed)
  d$y <- d$y * sc
  r <- fit_w(bf(y | mi() ~ x + ar(t, g, cov = TRUE)), family = lognormal(),
             data = d)
  dw <- d; dw$w <- ifelse(is.na(d$y), 0, 1); dw$y[is.na(dw$y)] <- sc
  rw <- fit_w(bf(y | weights(w) ~ x + ar(t, g, cov = TRUE)),
              family = lognormal(), data = dw)
  cat(sprintf("LN scale %g seed %d mi: %s | w0: %s\n", sc, seed, lat(r),
              lat(rw)))
}
gen <- function(seed, G = 20L, K = 10L, sd = 0.8) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  d$e <- as.vector(replicate(G, arima.sim(list(ar = 0.5), K, sd = sd)))
  d
}
for (seed in 1:3) {
  d <- gen(600 + seed)
  # gamma, shape 3 and 300
  d$yg <- rgamma(nrow(d), 3, rate = 3 / exp(0.5 + 0.3 * d$x + d$e))
  d$yg2 <- rgamma(nrow(d), 300, rate = 300 / exp(0.5 + 0.3 * d$x + d$e))
  # student with gaussian noise: nu large
  d$yt <- 1 + 0.3 * d$x + d$e + rnorm(nrow(d), 0, 0.5)
  d$yt2 <- (1 + 0.3 * d$x + d$e + rt(nrow(d), 3) * 0.5) * 1e6
  # weibull
  d$yw <- rweibull(nrow(d), 2, exp(0.5 + 0.3 * d$x + d$e))
  for (v in list(list("GAMMA3", bf(yg ~ x + ar(t, g, cov = TRUE)),
                      Gamma(link = "log")),
                 list("GAMMA300", bf(yg2 ~ x + ar(t, g, cov = TRUE)),
                      Gamma(link = "log")),
                 list("STUDgauss", bf(yt ~ x + ar(t, g, cov = TRUE)),
                      student()),
                 list("STUDt3e6", bf(yt2 ~ x + ar(t, g, cov = TRUE)),
                      student()),
                 list("WEIB", bf(yw ~ x + ar(t, g, cov = TRUE)), weibull()))) {
    r <- fit_w(v[[2]], family = v[[3]], data = d)
    cat(sprintf("%s seed %d: %s\n", v[[1]], seed, lat(r)))
  }
}
# few rows: 4 groups of 5
for (seed in 1:4) {
  d <- gen(700 + seed, G = 4L, K = 5L)
  d$yp <- rpois(nrow(d), exp(0.3 * d$x + d$e))
  d$yb <- rbinom(nrow(d), 1, plogis(0.3 * d$x + d$e))
  rp <- fit_w(bf(yp ~ x + ar(t, g, cov = TRUE)), family = poisson(), data = d)
  rb <- fit_w(bf(yb ~ x + ar(t, g, cov = TRUE)), family = bernoulli(),
              data = d)
  cat(sprintf("FEWPOIS %d: %s\nFEWBERN %d: %s\n", seed, lat(rp), seed, lat(rb)))
}
# poisson, identity link, mostly zeros: rows near mu = 0 leave the
# support at eta - 2h, and the gate skips non-finite rows
for (seed in 1:4) {
  d <- gen(800 + seed, sd = 0.3)
  mu <- pmax(0.05 + 0.02 * d$x + 0.1 * d$e, 0.01)
  d$yi <- rpois(nrow(d), mu)
  r <- fit_w(bf(yi ~ x + ar(t, g, cov = TRUE)), family = poisson("identity"),
             data = d)
  cat(sprintf("POISID %d: %s\n", seed, lat(r)))
}
cat("END\n")
