# Lane ci691: how far each xbeta precise-response fit sits from the SE
# check's flat threshold: the smallest over the largest eigenvalue of
# the unit-diagonal AD Hessian at the estimate, against se_flat_tol
# (1e-9), over the start perturbations of dev/ci691-xbeta-perturb.R.
#   Rscript dev/ci691-xbeta-ratio.R <lib> <phi> <kmax>
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r7"
phi <- if (length(args) > 1) as.numeric(args[2]) else 1e4
kmax <- if (length(args) > 2) as.integer(args[3]) else 24L
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r7",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "phi", phi, "\n")
set.seed(4200)
z <- stats::rbeta(1000, 0.047 * phi, 0.953 * phi)
d <- data.frame(y = pmin(pmax(1.1 * z - 0.05, 0), 1))
u <- suppressWarnings(frm(y ~ 1, family = xbeta(), data = d,
                          dry_run = "objective"))
p0 <- u$obj$par
r <- numeric(0)
for (k in 0:kmax) {
  st <- p0
  if (k > 0) {
    set.seed(k)
    st <- p0 * (1 + sample(c(-4:-1, 1:4), length(p0), TRUE) *
                  .Machine$double.eps)
  }
  fit <- suppressWarnings(frm(y ~ 1, family = xbeta(), data = d,
    start = list(beta = unname(st[1]), betad = unname(st[2:3]))))
  H <- fit$obj$he(fit$opt$par)
  H <- (H + t(H)) / 2
  D <- sqrt(abs(diag(H)))
  ev <- eigen(H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  r <- c(r, min(ev) / max(abs(ev)))
}
cat("fits", length(r), " min/max eigenvalue ratio: smallest",
    format(min(r), digits = 4), " largest", format(max(r), digits = 4),
    " se_flat_tol", frmtmb:::se_flat_tol, "\n")
