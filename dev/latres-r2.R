# R2: is M3's message true, on the four shapes of the review. Reads the build under
# test; the designs are dev/review-latres-m35b.R, -m3c.R and -m3d.R.
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
options(warn = 1, width = 120)
capt <- function(expr) {
  m <- w <- character(0)
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    message = function(c) { m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage") },
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  list(v = v, m = m, w = w)
}
txt <- function(r) paste(c(r$m, r$w), collapse = " | ")
cat("== ar at -1, alternating group effects, seeds 1 to 40\n")
tab <- list()
for (s in 1:40) {
  set.seed(s)
  d <- expand.grid(t = 1:6, g = factor(1:10))
  d$x <- rnorm(60)
  re <- as.vector(replicate(10, rnorm(1, 0, .6) * rep(c(1, -1), length.out = 6)))
  d$y <- rpois(60, exp(.5 + .3 * d$x + re))
  r <- capt(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
  f <- r$v
  th <- f$estimates$theta
  ar <- th[2] / sqrt(1 + th[2]^2)
  tab[[s]] <- data.frame(seed = s, ar = ar, cosy = grepl("read as cosy", txt(r)),
                         alt = grepl("alternates", txt(r)),
                         ll = as.numeric(logLik(f)))
}
tab <- do.call(rbind, tab)
cat("ar < -0.99:", sum(tab$ar < -.99), " says cosy:", sum(tab$cosy), " of those at -1:",
    sum(tab$cosy & tab$ar < 0), " says alternates:", sum(tab$alt), "\n")
cat("\n== a strong AR on a long series\n")
for (sd in 1:3) for (n in c(2000, 8000)) {
  set.seed(sd)
  e <- numeric(n)
  e[1] <- rnorm(1, 0, .05 / sqrt(1 - .9995^2))
  for (i in 2:n) e[i] <- .9995 * e[i - 1] + rnorm(1, 0, .05)
  d <- data.frame(t = seq_len(n), g = factor(1), x = rnorm(n))
  d$y <- rpois(n, exp(1 + .3 * d$x + e))
  r <- capt(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
  f <- r$v
  tt <- f$estimates$theta
  sm <- summary(f)$cor_pars
  cat(sprintf("seed %d n %d: ar %.5f sderr %.4f Est.Error ar %s sderr %s cosy %s boundary %s\n", sd, n,
              tt[2] / sqrt(1 + tt[2]^2), exp(tt[1]), format(sm[1, 2], digits = 3), format(sm[2, 2], digits = 3),
              grepl("read as cosy", txt(r)), grepl("Boundary", txt(r))))
}
cat("\n== skew_normal in three units, seed 7\n")
for (u in c(1, 1e-2, 1e-4)) {
  set.seed(7)
  dd <- expand.grid(t = 1:8, g = factor(1:25))
  dd <- dd[order(dd$g, dd$t), ]
  dd$x <- rnorm(nrow(dd))
  e <- as.vector(replicate(25, arima.sim(list(ar = .6), 8, sd = .5)))
  dd$y <- u * (.3 * dd$x + e + rnorm(nrow(dd), 0, .3))
  r <- capt(frm(y ~ x + ar(t, g, cov = TRUE), family = skew_normal(), data = dd))
  f <- r$v
  if (inherits(f, "error")) {
    cat("units", u, "ERR", substr(conditionMessage(f), 1, 120), "\n")
    next
  }
  tt <- f$estimates$theta
  sm <- summary(f)$cor_pars
  cat(sprintf("units %g: sderr/units %.3f ar %.3f Est.Error sderr %s says no latent %s\n", u, exp(tt[1]) / u,
              tt[2] / sqrt(1 + tt[2]^2), format(sm[2, 2], digits = 3), grepl("no latent residual", txt(r))))
}
cat("\n== arma on a random-intercept DGP, seeds 1 to 25\n")
for (s in 1:25) {
  set.seed(s)
  d <- expand.grid(t = 1:6, g = factor(1:10))
  d$x <- rnorm(60)
  d$y <- rpois(60, exp(.5 + .3 * d$x + rep(rnorm(10, 0, .5), each = 6)))
  r <- capt(frm(y ~ x + arma(t, g, cov = TRUE), family = poisson(), data = d))
  f <- r$v
  bk <- ns$frame_latres(f$frame)[[1]]
  lim <- ns$latres_block_limit(f, bk)
  sm <- summary(f)$cor_pars
  nan_ma <- !is.finite(sm["ma[1]", 2])
  named <- grepl("ma[1]", txt(r), fixed = TRUE)
  if (any(lim$at) || nan_ma) {
    cat(sprintf("seed %2d at %s ar %.4f ma %.3f sderr %.2g Est.Error ar %s ma %s sderr %s ma NaN or NA %s named %s cosy %s\n", s,
                paste(as.integer(lim$at), collapse = ""), sm["ar[1]", 1], sm["ma[1]", 1], sm["sderr", 1],
                format(sm["ar[1]", 2], digits = 2), format(sm["ma[1]", 2], digits = 2), format(sm["sderr", 2], digits = 2),
                nan_ma, named, grepl("read as cosy", txt(r))))
  }
}
cat("\n== sderr at 0, seed 5 of the review: the message and diagnose()\n")
set.seed(5)
d <- expand.grid(t = 1:6, g = factor(1:10))
d$x <- rnorm(60)
d$y <- rpois(60, exp(.5 + .3 * d$x + rep(rnorm(10, 0, .5), each = 6)))
r <- capt(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
cat(r$m, sep = "\n")
print(summary(r$v)$cor_pars)
print(diagnose(r$v, quiet = TRUE)$singular)
