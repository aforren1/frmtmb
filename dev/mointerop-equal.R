## Fits with three or more mo() categories: the lane must give what
## rellib-r14 gives. Prints every estimate and the log-likelihood to 17
## significant digits, and fitted() with standard errors; the two arms'
## outputs are compared by dev/mointerop-equal-cmp.R.
##   Rscript dev/mointerop-equal.R <lane|base> > dev/mointerop-out/equal-<arm>.txt
arm <- commandArgs(TRUE)[1]
base <- c("C:/Users/adf44/source/r/rellib-r14",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mointerop-lib", base)
} else base)
suppressMessages(library(frmtmb))
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$z <- rnorm(n)
d$g <- factor(sample(1:8, n, replace = TRUE))
d$o <- factor(c("a", "b", "c", "d")[d$k], ordered = TRUE)
d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + 0.3 * d$z +
  rnorm(8)[d$g] * 0.4 + rnorm(n)
d$kg <- c(1L, 2L, 5L)[(d$k %% 3L) + 1L]
d$cnt <- rpois(n, exp(0.2 + 0.3 * d$k))
models <- list(
  mok = list(y ~ mo(k), gaussian()),
  mokz = list(y ~ mo(k) + z, gaussian()),
  mokxz = list(y ~ mo(k) * z, gaussian()),
  mog = list(y ~ mo(k) + (mo(k) | g), gaussian()),
  moo = list(y ~ mo(o) + z, gaussian()),
  gap = list(y ~ mo(kg), gaussian()),
  sig = list(bf(y ~ z, sigma ~ mo(k)), gaussian()),
  pois = list(cnt ~ mo(k), poisson()))
for (nm in names(models)) {
  m <- models[[nm]]
  fit <- q(frm(m[[1]], family = m[[2]], data = d))
  cat("MODEL", nm, "\n")
  cat("logLik", format(as.numeric(logLik(fit)), digits = 17), "\n")
  est <- unlist(fit$estimates)
  for (i in seq_along(est)) {
    cat("est", names(est)[i], format(est[[i]], digits = 17), "\n")
  }
  fe <- q(fixef(fit))
  for (i in seq_len(nrow(fe))) {
    cat("fixef", rownames(fe)[i], format(fe[i, 1:2], digits = 17), "\n")
  }
  nd <- d[1:6, ]
  ft <- q(fitted(fit, newdata = nd))
  cat("fitted", format(c(ft[, 1:2]), digits = 17), "\n")
}
