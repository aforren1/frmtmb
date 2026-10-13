# Lane ci691: optimality certificates for the xbeta precise-response
# fits that do not need a standard error, over the start perturbations
# of dev/ci691-xbeta-perturb.R. Per fit: code, logLik, SEs kept, the
# test's max |g| * se (central differences of the reference, se from
# vcov), the Newton decrement g' H^-1 g / 2 of the reference with the
# fit's AD Hessian, and the log-likelihood the reference's own nlminb
# gains when started at the fit's estimate.
#   Rscript dev/ci691-xbeta-cert.R <lib> <phi> <kmax>
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r7"
phi <- if (length(args) > 1) as.numeric(args[2]) else 2e4
kmax <- if (length(args) > 2) as.integer(args[3]) else 24L
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r7",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "phi", phi, "\n")
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
set.seed(4200)
z <- stats::rbeta(1000, 0.047 * phi, 0.953 * phi)
d <- data.frame(y = pmin(pmax(1.1 * z - 0.05, 0), 1))
f <- function(p) sum(ref_xbeta(d$y, plogis(p[1]), exp(p[2]), exp(p[3])))
u <- suppressWarnings(frm(y ~ 1, family = xbeta(), data = d,
                          dry_run = "objective"))
p0 <- u$obj$par
for (k in 0:kmax) {
  st <- p0
  if (k > 0) {
    set.seed(k)
    st <- p0 * (1 + sample(c(-4:-1, 1:4), length(p0), TRUE) *
                  .Machine$double.eps)
  }
  fit <- suppressWarnings(frm(y ~ 1, family = xbeta(), data = d,
    start = list(beta = unname(st[1]), betad = unname(st[2:3]))))
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  est <- c(qlogis(dp$mu[1]), log(dp$phi[1]), log(dp$kappa[1]))
  g <- vapply(1:3, function(j) {
    h <- replace(numeric(3), j, 1e-5)
    (f(est + h) - f(est - h)) / 2e-5
  }, 0)
  se <- sqrt(diag(vcov(fit, full = TRUE)))
  H <- fit$obj$he(fit$opt$par)
  nd <- drop(t(g) %*% solve(H, g)) / 2
  o <- nlminb(est, function(p) -f(p),
              control = list(eval.max = 2000, iter.max = 2000))
  gain <- -o$objective - f(est)
  cat(sprintf(paste("k %2d code %d ll %.6f se_kept %-5s g*se %8.3g",
                    "newton %9.3g ref_gain %9.3g ref_code %d\n"),
              k, fit$opt$convergence, as.numeric(logLik(fit)),
              all(is.finite(se)), max(abs(g) * se), nd, gain,
              o$convergence))
}
