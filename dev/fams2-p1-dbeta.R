# Punch 1, B1: the interior log density of xbeta with its argument on the
# tape. RTMB::dbeta() against two hand-written forms (lbeta from RTMB,
# and the package's lbeta_ad(), which xbeta uses), value and gradient in
# (z, a, b), against a 200-bit Rmpfr reference, phi from 1e0 to 1e6.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(Rmpfr)})
hand <- function(z, a, b) {
  (a - 1) * log(z) + (b - 1) * log1p(-z) - RTMB::lbeta(a, b)
}
cands <- list(
  rtmb = function(p) RTMB::dbeta(p[1], p[2], p[3], log = TRUE),
  hand = function(p) hand(p[1], p[2], p[3]),
  pkg = function(p) {
    (p[2] - 1) * log(p[1]) + (p[3] - 1) * log1p(-p[1]) -
      frmtmb:::lbeta_ad(p[2], p[3])
  })
tapes <- lapply(cands, function(f) {
  F <- MakeTape(f, c(0.3, 2, 3)); list(F = F, J = F$jacfun())
})
ref <- function(z, a, b) {
  Z <- mpfr(z, 200); A <- mpfr(a, 200); B <- mpfr(b, 200)
  v <- (A - 1) * log(Z) + (B - 1) * log1p(-Z) - (lgamma(A) + lgamma(B) - lgamma(A + B))
  g <- c((A - 1) / Z - (B - 1) / (1 - Z),
         log(Z) - digamma(A) + digamma(A + B),
         log1p(-Z) - digamma(B) + digamma(A + B))
  list(v = as.numeric(v), g = as.numeric(g))
}
rows <- list()
for (lphi in 0:6) for (mu in c(1e-6, 0.01, 0.3, 0.5, 0.9)) {
  phi <- 10^lphi; a <- mu * phi; b <- (1 - mu) * phi
  sd <- sqrt(mu * (1 - mu) / (phi + 1))
  for (z in unique(pmin(pmax(c(mu, mu - 2 * sd, mu + 2 * sd, 1e-6, 1 - 1e-6),
                             1e-9), 1 - 1e-9))) {
    r <- ref(z, a, b)
    for (nm in names(tapes)) {
      tp <- tapes[[nm]]
      v <- tp$F(c(z, a, b)); g <- tp$J(c(z, a, b))
      rows[[length(rows) + 1L]] <- data.frame(
        phi = phi, mu = mu, z = z, form = nm, ref = r$v,
        verr = abs(v - r$v) / max(1, abs(r$v)),
        gfin = all(is.finite(g)),
        gerr = if (all(is.finite(g))) max(abs(g - r$g) / pmax(1, abs(r$g))) else NA)
    }
  }
}
res <- do.call(rbind, rows)
res <- res[is.finite(res$ref), ]
cat("points per form:", nrow(res) / 2, "\n")
for (nm in names(tapes)) {
  s <- res[res$form == nm, ]
  cat(sprintf("%-5s gradient finite %d of %d\n", nm, sum(s$gfin), nrow(s)))
  print(aggregate(cbind(verr, gerr) ~ phi, data = s,
                  FUN = function(x) signif(max(x), 2), na.action = na.pass))
}
