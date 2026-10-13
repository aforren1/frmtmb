# Lane ci691: the assertions the revised test-xbeta-zibb-hurdle-cum.R
# ("xbeta fits a precise response, with its ends") makes at phi 2e4,
# asked of fits from starts a few ulps apart, so that both outcomes of
# the SE check (kept and lost) are seen to pass. Per fit: SEs lost, the
# shown-or-named check, and the Newton gain against grad_tol.
#   Rscript dev/ci691-xbeta-testbody.R <lib> <phi> <kmax>
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
ref <- function(y, mu, phi, kappa) sum(ref_xbeta(y, mu, phi, kappa))
set.seed(4200)
z <- stats::rbeta(1000, 0.047 * phi, 0.953 * phi)
d <- data.frame(y = pmin(pmax(1.1 * z - 0.05, 0), 1))
u <- suppressWarnings(frm(y ~ 1, family = xbeta(), data = d,
                          dry_run = "objective"))
p0 <- u$obj$par
ok <- 0L
nlost <- 0L
gains <- numeric(0)
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
  ll <- as.numeric(stats::logLik(fit))
  r <- ref(d$y, dp$mu[1], dp$phi[1], dp$kappa[1])
  a1 <- abs(ll - r) < 1e-8 * abs(r)
  se <- sqrt(diag(vcov(fit, full = TRUE)))
  lost <- as.character(names(frmtmb:::sdr_of(fit)$se_lost))
  a2 <- setequal(names(se)[!is.finite(se)], lost)
  est <- c(stats::qlogis(dp$mu[1]), log(dp$phi[1]), log(dp$kappa[1]))
  g <- vapply(1:3, function(j) {
    h <- replace(numeric(3), j, 1e-5)
    f <- function(p) ref(d$y, stats::plogis(p[1]), exp(p[2]), exp(p[3]))
    (f(est + h) - f(est - h)) / 2e-5
  }, 0)
  H <- fit$obj$he(fit$opt$par)
  gain <- abs(sum(g * solve(H, g))) / 2
  a3 <- gain < fit$control$grad_tol
  ok <- ok + (a1 && a2 && a3)
  nlost <- nlost + (length(lost) > 0)
  gains <- c(gains, gain)
  cat(sprintf(paste("k %2d code %d lost %d [%s] ll_ok %s",
                    "shown_or_named %s gain %.3g\n"),
              k, fit$opt$convergence, length(lost),
              paste(lost, collapse = ","), a1, a2, gain))
}
cat(sprintf("fits %d  all assertions hold %d  SEs lost %d  max gain %.3g\n",
            kmax + 1L, ok, nlost, max(gains)))
