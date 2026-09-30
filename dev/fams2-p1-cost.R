# SUPERSEDED in punch 2: written against the punch-1
# log_ibeta_cf(x, a, b, N, lx, l1mx), whose signature punch 2 changed,
# so it no longer runs on the lane build. Its output is kept. The
# punch-2 cost is dev/fams2-p2-rev2-mean.txt, section 4.
# Punch 1, B2: what the new log incomplete beta costs inside a tape.
# Three versions of sum_i log I_{x_i}(a_i, b_i) over 500 boundary rows,
# each taped once in (x, log a, log b): the round-1 function (hard switch,
# fraction only; copied here verbatim), the new one, and RTMB::pbeta() as
# the floor. Reported: the time of one gradient sweep,
# the minimum over 7 rounds of blocks longer than 1.2 s, interleaved,
# with a control (the round-1 function taped twice) that must read 1.0.
# Seed 5.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
old <- function(x, a, b, N = 50L) {
  m <- (a + 1) / (a + b + 2)
  t <- m - x
  xd <- x - 0.5 * ((x - m) + abs(x - m))
  xc <- x + 0.5 * ((m - x) + abs(m - x))
  w <- (t + abs(t)) / (2 * abs(t) + 1e-300)
  ld <- frmtmb:::log_ibeta_cf(xd, a, b, N)
  lc <- frmtmb:::log_ibeta_cf(1 - xc, b, a, N)
  cap <- log1p(-2^-53)
  lc <- 0.5 * (lc + cap - abs(cap - lc))
  w * ld + (1 - w) * log1p(-exp(lc))
}
set.seed(5)
n <- 500
a <- exp(runif(n, log(0.5), log(5e3)))
b <- exp(runif(n, log(0.5), log(5e3)))
x <- pmin(runif(n, 0.01, 0.49), 0.49)
p0 <- c(x, log(a), log(b))
mk <- function(f) {
  MakeTape(function(p) {
    sum(f(p[seq_len(n)], exp(p[n + seq_len(n)]), exp(p[2 * n + seq_len(n)])))
  }, p0)
}
tapes <- list(control = mk(old), old = mk(old),
              new = mk(frmtmb:::log_ibeta_half),
              pbeta = mk(function(x, a, b) log(RTMB::pbeta(x, a, b))))
J <- lapply(tapes, function(F) F$jacfun())
timeit <- function(g) {
  k <- 1L
  repeat {
    # a fresh point each call, so no evaluation can be served from a cache
    t <- system.time(for (i in seq_len(k)) {
      g(p0 * (1 + 1e-9 * (i %% 7)))
    })[["elapsed"]]
    if (t > 1.2) return(t / k)
    k <- k * 2L
  }
}
res <- matrix(NA, 7, length(J), dimnames = list(NULL, names(J)))
for (r in 1:7) for (nm in names(J)) res[r, nm] <- timeit(J[[nm]])
best <- apply(res, 2, min)
cat("rows", n, "\n")
cat(sprintf("gradient sweep, ms (min of 7): %s\n",
            paste(sprintf("%s %.3f", names(best), 1e3 * best), collapse = "  ")))
cat(sprintf("ratios: control/old %.3f  new/old %.3f  new/pbeta %.3f\n",
            best[["control"]] / best[["old"]], best[["new"]] / best[["old"]],
            best[["new"]] / best[["pbeta"]]))
