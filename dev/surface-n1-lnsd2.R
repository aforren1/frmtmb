# Reviewer re-check, B2: where does logitnormal_sd()'s integrate() path
# fail, and the p = 1e-30, s = 10 discrepancy against a careful
# reference (fine pieces over the logistic's transition).
source("dev/surface-env.R"); surface_env("lane")
suppressPackageStartupMessages(library(frmtmb))
lnsd <- frmtmb:::logitnormal_sd
for (s in c(1e4, 3e4, 1e5, 3e5, 1e6, 1e7)) for (p in c(0.5, 1e-3, 1e-30)) {
  r <- tryCatch(format(lnsd(p, 1 - p, s), digits = 6),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("s %-6g p %-6g -> %s\n", s, p, r))
}
# careful reference for mu = qlogis(1e-30), s = 10: in z, pieces of
# width 0.5 over [mu - 14 s, mu + 14 s]
careful <- function(mu, s) {
  br <- seq(mu - 14 * s, mu + 14 * s, by = 0.5)
  I <- function(g) sum(vapply(seq_len(length(br) - 1L), function(j)
    stats::integrate(function(z) g(z) * stats::dnorm(z, mu, s), br[j],
                     br[j + 1], rel.tol = 1e-12, abs.tol = 0)$value, 0))
  m <- I(stats::plogis); sqrt(I(function(z) (stats::plogis(z) - m)^2))
}
for (cc in list(c(1e-30, 10), c(1e-30, 5), c(1e-20, 8), c(1e-12, 6),
                c(1e-8, 10))) {
  mu <- qlogis(cc[1]); s <- cc[2]
  a <- lnsd(cc[1], 1 - cc[1], s); b <- careful(mu, s)
  cat(sprintf("p %-6g s %-3g: frmtmb %-12.6g careful %-12.6g rel %.3g\n",
              cc[1], s, a, b, a / b - 1))
}
