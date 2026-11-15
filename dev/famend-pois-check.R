# Lane famend, 3580: poisson_cdf_log() with the Horner upper sum against
# the 0.82.0 form (git show HEAD:R/families-cdf.R, saved by the caller
# to dev/famend-out/families-cdf-head.R) and stats::ppois(), on a grid
# of counts 1 to 30 and means 1e-3 to 160: values, and the first three
# derivatives in log(mu) on the tape.
#   Rscript dev/famend-pois-check.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
src <- readLines("dev/famend-out/families-cdf-head.R")
i0 <- grep("^poisson_cdf_log <- function", src)
i1 <- i0 - 1L + grep("^}", src[i0:length(src)])[1]
old <- eval(parse(text = src[i0:i1]))
environment(old) <- ns
new <- ns$poisson_cdf_log
g <- expand.grid(q = 1:30, m = 10^seq(-3, log10(160), length.out = 80))
on <- new(g$q, g$m)
oo <- old(g$q, g$m)
rF <- stats::ppois(g$q, g$m, log.p = TRUE)
rS <- stats::ppois(g$q, g$m, lower.tail = FALSE, log.p = TRUE)
relerr <- function(x, r) {
  ok <- is.finite(r) & r != 0
  max(abs(x[ok] - r[ok]) / abs(r[ok]))
}
cat(sprintf("VALUE grid %d: new vs old max rel lF %.3g lS %.3g\n", nrow(g),
            relerr(on$lF, oo$lF), relerr(on$lS, oo$lS)))
cat(sprintf("VALUE new vs stats::ppois max rel lF %.3g lS %.3g; old %.3g %.3g\n",
            relerr(on$lF, rF), relerr(on$lS, rS), relerr(oo$lF, rF),
            relerr(oo$lS, rS)))
# where the Horner arm's mean is held: the survivor there
cat(sprintf("HELD max over q of P(Y <= q) at mu = q + 41: %.3g\n",
            max(stats::ppois(1:30, (1:30) + 41))))
bad <- 0
dif <- 0
for (i in seq_len(nrow(g))) {
  f <- function(fn) function(p) {
    l <- fn(g$q[i], exp(p)); l$lF + 2 * l$lS
  }
  x0 <- log(g$m[i])
  t1 <- MakeTape(f(new), x0)
  t0 <- MakeTape(f(old), x0)
  d1 <- c(t1$jacobian(x0), t1$jacfun()$jacobian(x0),
          t1$jacfun()$jacfun()$jacobian(x0))
  d0 <- c(t0$jacobian(x0), t0$jacfun()$jacobian(x0),
          t0$jacfun()$jacfun()$jacobian(x0))
  bad <- bad + any(!is.finite(d1))
  dif <- max(dif, abs(d1 - d0) / pmax(1, abs(d0)))
}
cat(sprintf(paste0("DERIV grid %d: rows with a non-finite first, second ",
                   "or third derivative %d; max rel |new - old| %.3g\n"),
            nrow(g), bad, dif))
