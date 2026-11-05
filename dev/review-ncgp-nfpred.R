# Item 2: are "the predictions right" on noise-free data, off the observed positions and
# in their standard errors? Default frm(), frm() started at the marginal's maximum, the
# marginal itself. Also which rule NaNs the intercept's SE at y * 1000.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
krig <- function(s, ux, y, b, x0) {
  C <- exp(-s$D2 / (2 * b$rho^2)); diag(C) <- diag(C) + s$nug + b$r
  R <- chol(C)
  Ks <- exp(-outer(x0, ux, function(a, c) (a - c)^2) / (2 * b$rho^2))
  al <- backsolve(R, backsolve(R, y - b$beta, transpose = TRUE))
  w <- backsolve(R, t(Ks), transpose = TRUE)
  list(mean = as.vector(b$beta + Ks %*% al), se = b$sd * sqrt(pmax(1 - colSums(w^2), 0)))
}
go <- function(lab, x, y, truef, x0) {
  s <- cf_setup(y, x)
  d <- data.frame(y = y, x = x)
  r <- rev_fit(y ~ gp(x), d); f <- r$fit; g <- rev_gp_read(f)
  m <- cf_max(s, extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 30L)
  ux <- x[match(unique(x), x)]
  ym <- rowsum_fast(y, s$J, s$m) / s$N
  kc <- krig(s, ux, ym, m, x0)
  pr <- suppressWarnings(fitted(f, newdata = data.frame(x = x0)))
  r2 <- rev_fit(y ~ gp(x), d, start = list(theta = c(log(m$sd), log(m$rho)), betad = log(m$sigma), beta = m$beta[1]))
  f2 <- r2$fit; g2 <- rev_gp_read(f2)
  pr2 <- suppressWarnings(fitted(f2, newdata = data.frame(x = x0)))
  tr <- truef(x0)
  cat(sprintf("RNF %s n %d | default: ll %.3f (cf %.3f, exact gap at point %.3f) sd %.3g rho %.3g sigma %.3g conv %d w %d | from cf: ll %.3f sd %.3g rho %.3g conv %d w %d\n",
              lab, length(y), g$ll, m$ll, cf_ll_mpfr(s, g$beta, g$sigma, g$sd, g$rho) - cf_ll_mpfr(s, m$beta, m$sigma, m$sd, m$rho), g$sd, g$rho, g$sigma, g$conv, length(r$warn), g2$ll, g2$sd, g2$rho, g2$conv, length(r2$warn)))
  cat(sprintf("RNF   at new x, max |pred - truth|: default %.2e from-cf %.2e marginal %.2e | se(default)/se(marginal) %s | se(from cf)/se(marginal) %s\n",
              max(abs(pr[, 1] - tr)), max(abs(pr2[, 1] - tr)), max(abs(kc$mean - tr)),
              paste(sprintf("%.3g", pr[, 2] / kc$se), collapse = " "), paste(sprintf("%.3g", pr2[, 2] / kc$se), collapse = " ")))
  flush.console()
}
set.seed(20261102)
n <- 60L
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
x0 <- c(0.507, 1.313, 2.223, 3.141, 4.444, 5.03, 5.9)
go("lane noise-free y", x, sin(2 * x), function(z) sin(2 * z), x0)
go("lane noise-free y*1000", x, 1000 * sin(2 * x), function(z) 1000 * sin(2 * z), x0)
go("lane noise-free y/1000", x, sin(2 * x) / 1000, function(z) sin(2 * z) / 1000, x0)
go("lane noise-free x*1000", x * 1000, sin(2 * x), function(z) sin(2 * z / 1000), x0 * 1000)
go("lane s=1e-8", x, sin(2 * x) + 1e-8 * e0, function(z) sin(2 * z), x0)
set.seed(21); x <- sort(runif(30, 0, 10))
go("n30 sin(x)", x, sin(x), sin, c(1.1, 3.3, 5.5, 7.7, 9.5))
cat("== the SE rule\n")
set.seed(14)
x <- sort(runif(40, 0, 8)); y <- sin(x) + 1e-4 * rnorm(40)
for (c0 in c(1, 300, 1000)) {
  f <- suppressWarnings(frm(y ~ gp(x), data.frame(y = y * c0, x = x)))
  sdr <- f$sdr %||% frmtmb:::sdr_of(f)
  H <- tryCatch(optimHess(f$opt$par, f$obj$fn, f$obj$gr), error = function(e) NULL)
  cat(sprintf("REQS c %g conv %d se_lost %s cov.fixed diag %s units %s | Hessian(gr) diag %s\n", c0, f$opt$convergence, paste(names(sdr$se_lost), sdr$se_lost, collapse = " "),
              paste(signif(diag(sdr$cov.fixed), 3), collapse = " "), paste(signif(f$par_units %||% NA, 3), collapse = " "), paste(signif(diag(H), 3), collapse = " ")))
  rawv <- sdr$cov_fixed_raw
  if (!is.null(rawv)) cat("REQS   raw diag", signif(diag(rawv), 4), "\n")
}
cat("DONE\n")
