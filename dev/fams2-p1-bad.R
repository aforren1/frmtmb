.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
source("C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-p1-cand.R")
set.seed(20260930); M <- 3000
a <- exp(runif(M, log(1e-3), log(1e7))); b <- exp(runif(M, log(1e-3), log(1e7)))
s <- a + b; m <- (a + 1) / (s + 2); sd <- sqrt(a * b / (s * s * (s + 1))); u <- runif(M)
x <- ifelse(u < 0.1, m, ifelse(u < 0.6, m + runif(M, -5, 5) * sd, runif(M, 0, 0.5)))
ok <- x > 0 & x < 0.5; x <- x[ok]; a <- a[ok]; b <- b[ok]
ref <- suppressWarnings(pbeta(x, a, b, log.p = TRUE)); keep <- is.finite(ref) & ref > -600
x <- x[keep]; a <- a[keep]; b <- b[keep]; ref <- ref[keep]
v <- mapply(cand, x, a, b)
bad <- which(!is.finite(v))
print(data.frame(x = x[bad], a = a[bad], b = b[bad], ref = ref[bad], v = v[bad]))
for (i in bad) {
  s <- a[i] + b[i]; m <- (a[i]+1)/(s+2); sd <- sqrt(a[i]*b[i]/(s*s*(s+1)))
  cat("k", (x[i]-m)/sd, " ld", frmtmb:::log_ibeta_cf(min(x[i], m+sd), a[i], b[i], 50L),
      " lc", frmtmb:::log_ibeta_cf(1 - max(x[i], m - sd), b[i], a[i], 50L),
      " pbeta", pbeta(x[i], a[i], b[i]), "\n")
}
dbg <- cand
body(dbg) <- as.call(c(as.list(body(cand))[-length(body(cand))], quote(print(c(wd = wd, ld = ld, lc = lc, lcf = lcf, k = k, ws = ws, wk = wk, wp = wp, x2 = x2, lp = lp)))))
dbg(x[bad[1]], a[bad[1]], b[bad[1]])
