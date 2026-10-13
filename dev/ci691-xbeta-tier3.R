# Lane ci691: why the SE check calls all three xbeta parameters flat at
# a converged phi 2e4 fit (start perturbation k of
# dev/ci691-xbeta-perturb.R), traced through se_tier3()'s quantities.
#   Rscript dev/ci691-xbeta-tier3.R <lib> <k>
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r7"
kk <- if (length(args) > 1) as.integer(args[2]) else 2L
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r7",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "k", kk, "\n")
phi <- 2e4
set.seed(4200)
z <- stats::rbeta(1000, 0.047 * phi, 0.953 * phi)
d <- data.frame(y = pmin(pmax(1.1 * z - 0.05, 0), 1))
u <- suppressWarnings(frm(y ~ 1, family = xbeta(), data = d,
                          dry_run = "objective"))
p0 <- u$obj$par
st <- p0
if (kk > 0) {
  set.seed(kk)
  st <- p0 * (1 + sample(c(-4:-1, 1:4), length(p0), TRUE) *
                .Machine$double.eps)
}
fit <- suppressWarnings(frm(y ~ 1, family = xbeta(), data = d,
  start = list(beta = unname(st[1]), betad = unname(st[2:3]))))
p <- fit$opt$par
cat("code", fit$opt$convergence, "par", format(p, digits = 10), "\n")
cat("gradient", format(fit$obj$gr(p), digits = 6), "\n")
H <- fit$obj$he(p)
H <- (H + t(H)) / 2
D <- sqrt(diag(H))
S <- H / outer(D, D)
e <- eigen(S, symmetric = TRUE)
cat("unit-diagonal eigenvalues", format(e$values, digits = 6), "\n")
cat("ratio min/max", format(min(e$values) / max(e$values), digits = 6),
    " se_flat_tol", frmtmb:::se_flat_tol, "\n")
cat("eigenvectors:\n"); print(e$vectors)
tol <- fit$control$grad_tol
cat("grad_tol", tol, "\n")
k <- 3
dir <- e$vectors[, k] / D
f <- function(x) fit$obj$fn(x)
for (s in c(1, 0.5, 0.25, 0.125)) {
  dd <- dir * sqrt(4 * tol / e$values[k]) * s
  c2 <- f(p + dd) + f(p - dd) - 2 * f(p)
  cat(sprintf("step %.3f  |d| %s  second difference %.6g  quad predicts %.6g\n",
              s, paste(format(dd, digits = 4), collapse = " "), c2,
              4 * tol * s^2))
}
cat("se_curvature_real:",
    frmtmb:::se_curvature_real(fit, p, 1:3, dir, e$values[k]), "\n")
ch <- frmtmb:::cov_from_hessian(fit, H, matrix(0, 3, 3), exact = TRUE)
cat("tier", ch$tier, "lost", paste(names(ch$lost), ch$lost), "\n")
cat("plain solve(H) se", format(sqrt(diag(solve(H))), digits = 6), "\n")
# the same direction read with the stats-function reference, so the
# objective's own rounding is out of the picture
ref_xbeta <- function(y, mu, phi, kappa) {
  a <- mu * phi
  b <- (1 - mu) * phi
  dd <- 1 + 2 * kappa
  ifelse(y <= 0, stats::pbeta(kappa / dd, a, b, log.p = TRUE),
         ifelse(y >= 1,
                stats::pbeta((1 + kappa) / dd, a, b, lower.tail = FALSE,
                             log.p = TRUE),
                stats::dbeta((y + kappa) / dd, a, b, log = TRUE) - log(dd)))
}
nll <- function(q) -sum(ref_xbeta(d$y, plogis(q[1]), exp(q[2]), exp(q[3])))
for (s in c(1, 0.5, 0.25)) {
  dd <- dir * sqrt(4 * tol / e$values[k]) * s
  c2 <- nll(p + dd) + nll(p - dd) - 2 * nll(p)
  cat(sprintf("reference: step %.3f second difference %.6g\n", s, c2))
}
for (t in c(-8, -4, -2, -1, 1, 2, 4, 8)) {
  dd <- dir * sqrt(4 * tol / e$values[k]) * t
  cat(sprintf("along the direction, %3g full steps: loss %.6g\n", t,
              f(p + dd) - f(p)))
}
