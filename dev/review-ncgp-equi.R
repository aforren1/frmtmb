# Item 2: equivariance. y * 1000, y / 1000, x * 1000, x / 1000 against the
# base design, in one process, at estimates, fitted values, SEs and kriging.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- commandArgs(TRUE)[1]
rev_libs(arm)
fitread <- function(y, x, xnew, ysc, xsc, fml = y ~ gp(x), ...) {
  d <- data.frame(y = y * ysc, x = x * xsc)
  r <- rev_fit(fml, d, ...)
  f <- r$fit
  if (inherits(f, "error")) return(list(err = conditionMessage(f)))
  g <- rev_gp_read(f)
  ft <- fitted(f)
  nd <- data.frame(x = xnew * xsc)
  pr <- fitted(f, newdata = nd)
  th <- sqrt(diag(vcov(f, full = TRUE)))
  list(ll = g$ll + length(y) * log(ysc), sd = g$sd / ysc, rho = g$rho / xsc,
       sigma = g$sigma / ysc, mu = g$beta / ysc, fit = ft[, 1] / ysc, se = ft[, 2] / ysc,
       pr = pr[, 1] / ysc, pse = pr[, 2] / ysc, conv = g$conv, nw = length(r$warn),
       ses = tryCatch(unname(sqrt(diag(as.matrix(vcov(f))))), error = function(e) NA))
}
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
for (case in list(list(s = 0.2, om = 1, seed = 11, n = 60), list(s = 1e-2, om = 2, seed = 12, n = 60),
                  list(s = 1e-3, om = 2, seed = 13, n = 80), list(s = 1e-4, om = 1, seed = 14, n = 40))) {
  set.seed(case$seed)
  x <- sort(runif(case$n, 0, 8))
  y <- sin(case$om * x) + case$s * rnorm(case$n)
  xn <- seq(0.3, 7.7, length.out = 7) + 0.0123
  b <- fitread(y, x, xn, 1, 1)
  for (sc in list(c(1000, 1), c(1e-3, 1), c(1, 1000), c(1, 1e-3), c(1e6, 1e-6))) {
    a <- fitread(y, x, xn, sc[1], sc[2])
    if (!is.null(a$err)) { cat("RE ERROR", a$err, "\n"); next }
    cat(sprintf(paste0("RE s %g om %g ysc %g xsc %g | conv %d/%d nw %d/%d | ll %.2e sd %.2e rho %.2e ",
                       "sigma %.2e mu %.2e fit %.2e fitse %.2e pred %.2e predse %.2e\n"),
                case$s, case$om, sc[1], sc[2], b$conv, a$conv, b$nw, a$nw,
                abs(a$ll - b$ll), rel(a$sd, b$sd), rel(a$rho, b$rho), rel(a$sigma, b$sigma),
                rel(a$mu, b$mu), rel(a$fit, b$fit), rel(a$se, b$se), rel(a$pr, b$pr), rel(a$pse, b$pse)))
    flush.console()
  }
}
cat("DONE\n")
