# Review lane osa: case weights under the "cdf" route.
#
#   Rscript dev/review-osa-weights.R [lib]
#
# A frequency weight of 2 (brms: the row's log-likelihood times 2) is the
# same fitted model as the row written twice. Each one-step interval
# [lo, hi] of the lane's residual is recovered from two calls with
# different seeds: pnorm(r) = hi - U (hi - lo) for each U.
# Exact one-step intervals by quadrature (1601 points over +-8 sd) at the
# fit's estimates, with the earlier rows' mass raised to their weight
# ("weighted conditioning", the fitted posterior of the group effect) or
# not ("unweighted conditioning").
args <- commandArgs(TRUE)
lib1 <- if (length(args) >= 1) args[1] else "C:/Users/adf44/source/r/wt-osa-lib"
.libPaths(c(if (lib1 != "none") lib1, "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-weights\n")

steps <- function(fit) {
  n <- nrow(fit$frame$data %||% fit$data)
  ra <- residuals(fit, type = "osa", seed = 123)[, "Estimate"]
  rb <- residuals(fit, type = "osa", seed = 456)[, "Estimate"]
  n <- length(ra)
  Ua <- local({ set.seed(123); runif(n) })
  Ub <- local({ set.seed(456); runif(n) })
  p <- (pnorm(ra) - pnorm(rb)) / (Ub - Ua)
  hi <- pnorm(ra) + Ua * p
  list(lo = hi - p, hi = hi)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

G <- 20L
K <- 10L
set.seed(301)
d <- data.frame(g = factor(rep(seq_len(G), each = K)), t = rep(1:K, G),
                x = rnorm(G * K))
d$y <- rpois(nrow(d), exp(0.2 + 0.5 * d$x + rnorm(G, 0, 1)[d$g]))
d$w <- ifelse(d$t <= 3, 2, 1)
# the same data with each weight-2 row written twice, the copy right after
dd <- d[rep(seq_len(nrow(d)), d$w), ]
dd$w <- 1
dd$orig <- rep(seq_len(nrow(d)), d$w)
dd$last <- !duplicated(dd$orig, fromLast = TRUE)

fw <- frm(y | weights(w) ~ x + (1 | g), family = poisson(), data = d)
fd <- frm(y ~ x + (1 | g), family = poisson(), data = dd)
cat(sprintf("estimates weighted %s | duplicated %s\n",
            paste(format(c(fw$estimates$beta, fw$estimates$theta), digits = 6),
                  collapse = " "),
            paste(format(c(fd$estimates$beta, fd$estimates$theta), digits = 6),
                  collapse = " ")))
sw <- steps(fw)
sd_ <- steps(fd)
# the original rows of the duplicated fit, each at its LAST copy
hi_d <- sd_$hi[dd$last]
lo_d <- sd_$lo[dd$last]
later <- d$t >= 4
cat(sprintf(paste("rows t >= 4 (weight 1, after the weighted rows):",
                  "max|qnorm(hi) weighted - duplicated| %.3g,",
                  "median %.3g\n"),
            max(abs(qnorm(sw$hi) - qnorm(hi_d))[later]),
            median(abs(qnorm(sw$hi) - qnorm(hi_d))[later])))

exact <- function(fit, d, wt) {
  eta <- fit$estimates$beta[[1]] + fit$estimates$beta[[2]] * d$x
  sdu <- exp(fit$estimates$theta[[1]])
  u <- seq(-8, 8, length.out = 1601) * sdu
  w0 <- dnorm(u, 0, sdu)
  lo <- hi <- numeric(nrow(d))
  for (gi in split(seq_len(nrow(d)), d$g)) {
    w <- w0 / sum(w0)
    for (i in gi) {
      m <- exp(eta[i] + u)
      p <- dpois(d$y[i], m)
      Fh <- ppois(d$y[i], m)
      hi[i] <- sum(w * Fh)
      lo[i] <- sum(w * (Fh - p))
      w <- w * p^(if (wt) d$w[i] else 1)
      w <- w / sum(w)
    }
  }
  list(lo = lo, hi = hi)
}
ew <- exact(fw, d, TRUE)
eu <- exact(fw, d, FALSE)
zq <- function(s) qnorm(s$hi)
cat(sprintf(paste("weighted fit, rows t >= 4: max|lane - exact weighted|",
                  "%.3g, max|lane - exact unweighted| %.3g,",
                  "max|exact weighted - exact unweighted| %.3g",
                  "(qnorm of the upper end)\n"),
            max(abs(zq(sw) - zq(ew))[later]), max(abs(zq(sw) - zq(eu))[later]),
            max(abs(zq(ew) - zq(eu))[later])))
U <- local({ set.seed(123); runif(nrow(d)) })
rw <- residuals(fw, type = "osa")[, "Estimate"]
zw <- qnorm(ew$hi - U * (ew$hi - ew$lo))
cat(sprintf(paste("weighted fit residual vs exact weighted-conditioning",
                  "residual, rows t >= 4: max %.3g median %.3g; rows over",
                  "0.05: %d of %d\n"),
            max(abs(rw - zw)[later]), median(abs(rw - zw)[later]),
            sum(abs(rw - zw)[later] > 0.05), sum(later)))

# gaussian: the default route (fullGaussian) against "cdf", weight 3
set.seed(302)
dg <- data.frame(x = rnorm(120))
dg$y <- 1 + 0.5 * dg$x + rnorm(120)
dg$w <- rep(c(1, 3), 60)
fg <- frm(y | weights(w) ~ x, family = gaussian(), data = dg)
r_fg <- residuals(fg, type = "osa")[, "Estimate"]
r_cdf <- residuals(fg, type = "osa", osa_method = "cdf")[, "Estimate"]
fs <- frm(y | weights(w) ~ x, family = student(), data = dg)
r_st <- residuals(fs, type = "osa")[, "Estimate"]
cat(sprintf(paste("gaussian weight 3 rows: mean r_fullGaussian / r_cdf %.4f",
                  "(sqrt(3) = %.4f); weight 1 rows: %.4f.",
                  "student() same data, sd of weight-3 rows %.3f vs",
                  "gaussian default %.3f\n"),
            mean(r_fg[dg$w == 3] / r_cdf[dg$w == 3]), sqrt(3),
            mean(r_fg[dg$w == 1] / r_cdf[dg$w == 1]),
            sd(r_st[dg$w == 3]), sd(r_fg[dg$w == 3])))
cat("END review-osa-weights\n")
