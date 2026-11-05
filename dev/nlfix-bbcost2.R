# R3: the shipped families, whole fits, on one library at a time. Arg 1 is the
# library. Alternate the two arms, since the machine is shared; read the minima.
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")), "\n")
mk <- function(n, ntr, phi = 8) {
  set.seed(77)
  d <- data.frame(x = rnorm(n), ntr = as.integer(ntr))
  p <- plogis(0.3 + 0.5 * d$x)
  d$y <- rbinom(n, d$ntr, rbeta(n, p * phi, (1 - p) * phi))
  d$z <- d$y * (runif(n) > 0.1)
  d$cens <- ifelse(runif(n) < 0.2 & d$y < d$ntr, "right", "none")
  d
}
tm <- function(f, k = 4L) {
  r <- NULL
  v <- numeric(k)
  for (i in seq_len(k)) {
    gc(FALSE)
    v[i] <- system.time(r <- suppressWarnings(suppressMessages(f())))[["elapsed"]]
  }
  c(min(v), as.numeric(logLik(r)), r$opt$convergence)
}
for (ntr in if (length(args) > 1L) numeric(0) else c(5, 10, 24, 30, 200)) {
  d <- mk(3e4, ntr)
  a <- tm(function() frm(bf(y | trials(ntr) ~ x), family = beta_binomial(), data = d))
  b <- tm(function() frm(bf(z | trials(ntr) ~ x), family = zero_inflated_beta_binomial(), data = d), 3L)
  cat(sprintf("trials %3d: bb %6.2f s ll %.6f code %d | zibb %6.2f s ll %.6f code %d\n",
              ntr, a[1], a[2], a[3], b[1], b[2], b[3]))
}
d <- mk(3000, 10)
a <- tm(function() frm(bf(y | trials(ntr) + cens(cens) ~ x), family = beta_binomial(), data = d), 3L)
cat(sprintf("censored, n = 3000, trials 10: %.2f s ll %.6f code %d\n", a[1], a[2], a[3]))
d <- mk(3000, 200)
a <- tm(function() frm(bf(y | trials(ntr) + cens(cens) ~ x), family = beta_binomial(), data = d), 3L)
cat(sprintf("censored, n = 3000, trials 200: %.2f s ll %.6f code %d\n", a[1], a[2], a[3]))
