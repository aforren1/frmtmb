.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-lib", "C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
sim_hnb <- function(seed, n = 500, shape = 1.3) {
  set.seed(seed); x <- rnorm(n); z <- rnorm(n)
  mu <- exp(0.6 + 0.4 * x); hu <- plogis(-0.4 + 0.7 * z)
  p0 <- stats::dnbinom(0, size = shape, mu = mu)
  yp <- pmax(stats::qnbinom(p0 + runif(n) * (1 - p0), size = shape, mu = mu), 1)
  data.frame(y = ifelse(runif(n) < hu, 0, yp), x = x, z = z)
}
d <- sim_hnb(5, n = 60)
r <- try(frm(bf(y | trunc(lb = 1) ~ x) + hurdle_negbinomial(), data = d))
r2 <- try(frm(bf(y | trunc(lb = 1) ~ x) + poisson(), data = d))
r3 <- try(frm(bf(y | cens(x > 1) ~ x) + hurdle_negbinomial(), data = d))
print(r3)
d2 <- d[d$y >= 1, ]
r4 <- try(frm(bf(y | trunc(lb = 1) ~ x) + hurdle_negbinomial(), data = d2, verbose = FALSE))
print(r4)
