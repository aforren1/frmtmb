# Fix round 2: the upper tail of the CDF, which dev/review-nlfix-d3.R did not print; and make
# dev/review-nlfix-c1.R section H run on this build.
p <- "dev/review-nlfix-c1.R"; l <- readLines(p); i <- grep("nl_restart_credible", l)
if (length(i)) { l[i] <- "  # nl_restart_credible() was removed in fix round 2"; writeLines(l, p) }
source("dev/review-nlfix-common.R")
ns <- asNamespace("frmtmb")
fp <- function(y, n, a, b) { phi <- a + b; lchoose(n, y) + (if (y > 0) sum(log(a + 0:(y - 1))) else 0) + (if (n - y > 0) sum(log(b + 0:(n - y - 1))) else 0) - sum(log(phi + 0:(n - 1))) }
worst <- 0; bad <- 0L
for (n in c(1, 5, 10, 24, 25, 60)) for (lp in c(-3, 0, 2, 15, 30, 39, 60, 100)) for (m in c(0.02, 0.5, 0.98)) {
  phi <- exp(lp); qs <- sort(unique(c(-1, 0, 1, floor(n / 2), n - 1, n, n + 1))); N <- length(qs)
  cd <- ns$beta_binomial_cdf_log(qs, rep(n, N), rep(m * phi, N), rep((1 - m) * phi, N))
  hs <- vapply(qs, function(q) { if (q < 0) return(0); if (q >= n) return(-Inf); log(sum(exp(sapply((q + 1):n, fp, n = n, a = m * phi, b = (1 - m) * phi)))) }, 0)
  v <- as.numeric(cd$lS)
  if (any(is.finite(hs) != is.finite(v))) { bad <- bad + 1L; if (bad <= 3) cat(" mismatch n", n, "lp", lp, "m", m, ":", v, "|", hs, "\n") }
  ok <- is.finite(hs) & is.finite(v); worst <- max(worst, abs(v[ok] - hs[ok]) / pmax(1, abs(hs[ok])))
}
cat("upper tail lS: largest relative difference", format(worst, digits = 3), " finite mismatches", bad, "\n")
cat("PID 7544:", system2("tasklist", c("/FI", shQuote("PID eq 7544"), "/NH"), stdout = TRUE)[1:2], sep = "\n")
