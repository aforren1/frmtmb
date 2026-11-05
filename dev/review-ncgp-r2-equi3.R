# Fix round 1, item 2: y and x scaled on fresh curves, not dev/ncgp-equiv.R's. Gap of each scale against scale 1 and
# against the closed form; NaN SEs. Seeds 62000 + n.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
fns <- list(tan1 = function(x) atan(5 * (x - 5)) + cos(9 * x) / 2, poly = function(x) (x - 2) * (x - 5) * (x - 8) / 10 + sin(4 * x), gauss4 = function(x) sin(7 * x) * exp(-x / 5) + 0.3 * x, sin6 = function(x) sin(6 * x))
tot <- 0; bad <- 0; nan <- 0
for (fn in names(fns)) for (sg in c(0.1, 1e-3, 1e-4)) for (n in c(40L, 90L)) {
  set.seed(62000 + n); x <- sort(runif(n, 0, 10)); e <- rnorm(n); y <- fns[[fn]](x) + sg * e
  nd <- c(0.33, 4.44, 8.88)
  rd <- function(yc, xc) {
    r <- rev_fit(y ~ gp(x), data.frame(y = y * yc, x = x * xc)); g <- rev_gp_read(r$fit)
    p <- suppressWarnings(fitted(r$fit, newdata = data.frame(x = nd * xc)))
    list(ll = g$ll + n * log(yc), sd = g$sd / yc, rho = g$rho / xc, sigma = g$sigma / yc, se = p[, 2] / yc, conv = g$conv, w = length(r$warn))
  }
  b <- rd(1, 1)
  m <- cf_max(cf_setup(y, x), ngrid = 22L)
  for (sc in list(c(1000, 1), c(1e-3, 1), c(1e6, 1), c(1, 1000), c(1, 1e-3))) {
    a <- rd(sc[1], sc[2])
    tot <- tot + 1; if (abs(a$ll - b$ll) > 1e-3) bad <- bad + 1; if (any(is.nan(a$se))) nan <- nan + 1
    cat(sprintf("RE3Q %s s=%g n %d y*%g x*%g | ll vs scale 1 %.2e, vs cf %.2e, scale 1 vs cf %.2e | sd %.1e rho %.1e sigma %.1e se %.1e | conv %d w %d\n", fn, sg, n, sc[1], sc[2], a$ll - b$ll, a$ll - m$ll, b$ll - m$ll,
                abs(a$sd / b$sd - 1), abs(a$rho / b$rho - 1), abs(a$sigma / b$sigma - 1), max(abs(a$se / b$se - 1)), a$conv, a$w))
    flush.console()
  }
}
cat(sprintf("RE3QSUM rows %d differ by more than 1e-3 nats %d NaN SE %d\n", tot, bad, nan))
cat("DONE\n")
