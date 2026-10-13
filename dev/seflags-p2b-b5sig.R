# Lane seflags, punch round 2b: the B5 test's sigma rows (seeds 3, 4, 10
# with (1 | g)) against poly(x, 8), with each fit's sigma, its
# log-likelihood against the reference's, and the ill-conditioning
# warning (which carries the Newton gain).
#   Rscript dev/seflags-p2b-b5sig.R <lane|m9>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev3-common.R")
for (s in c(3, 4, 10)) {
  set.seed(s)
  d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)),
                  z = stats::rnorm(200), w = factor(sample(c("a", "b"), 200, TRUE)))
  d$y <- sin(3 * d$x) + 0.3 * d$z + 0.5 * (d$w == "b") +
    stats::rnorm(20, 0, 0.3)[d$g] + stats::rnorm(200, 0, 0.2)
  for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
  r <- cap(frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + z + w + (1 | g), data = d))
  f <- r$v
  fr <- suppressWarnings(frm(y ~ poly(x, 8) + z + w + (1 | g), data = d))
  nd <- data.frame(x = 1.5, z = 0, w = "a", g = d$g[1])
  for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
  a <- suppressWarnings(fitted(f, newdata = nd, re_formula = NA, dpar = "sigma"))
  b <- suppressWarnings(fitted(fr, newdata = nd, re_formula = NA, dpar = "sigma"))
  cat(sprintf("seed %d code %d | sigma %.5f ref %.5f | SE %.5f ref %.5f ratio %.4f | logLik diff %.4f | %s\n",
              s, f$opt$convergence, a[1, 1], b[1, 1], a[1, 2], b[1, 2], a[1, 2] / b[1, 2],
              as.numeric(logLik(f)) - as.numeric(logLik(fr)),
              sub(".*(Newton step there raises the log-likelihood by [0-9.e-]+).*", "\1",
                  paste(r$w, collapse = " "))))
}
