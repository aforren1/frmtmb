source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages(library(frmtmb))
set.seed(20260917)
d <- data.frame(g = factor(rep(1:6, each = 10)), x = rnorm(60))
d$y <- d$x + rnorm(60)
fit <- frm(y ~ x + (1 | g), d)
alt <- fit
class(alt) <- c("revdollar", class(alt))
registerS3method("$", "revdollar", function(x, name) {
  .subset2(x, name, exact = FALSE)
})
per <- function(obj, n, rounds = 5) {
  best <- Inf
  for (r in seq_len(rounds)) {
    t0 <- proc.time()[["elapsed"]]
    for (i in seq_len(n)) obj$spec
    el <- proc.time()[["elapsed"]] - t0
    best <- min(best, el)
  }
  list(us = best / n * 1e6, secs = best, n = n)
}
a <- per(fit, 4e6); a2 <- per(fit, 4e6); b <- per(alt, 1e6)
cat(sprintf("plain    %.3f us (block %.2f s, n=%g)\n", a$us, a$secs, a$n))
cat(sprintf("plain2   %.3f us (block %.2f s, n=%g)\n", a2$us, a2$secs, a2$n))
cat(sprintf("method   %.3f us (block %.2f s, n=%g)\n", b$us, b$secs, b$n))
cat(sprintf("CONTROL plain2/plain = %.4f (must be about 1.0)\n",
            a2$us / a$us))
cat(sprintf("method - plain = %.3f us; at 751 reads = %.3f ms\n",
            b$us - a$us, 751 * (b$us - a$us) / 1000))
