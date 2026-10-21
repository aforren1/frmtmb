# Release reviewer, 0.73.0: peak memory of a single-block quadrature fit
# with a censored response, by family and size; the caller
# (dev/rel073-rev-capped.ps1) records the process's peak working set.
#   Rscript dev/rel073-rev-quadmem.R <lib> <family> <n> <quadrature TRUE|FALSE>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
fam <- a[2]; n <- as.integer(a[3]); quad <- as.logical(a[4])
cat("lib: frmtmb", format(packageVersion("frmtmb")), "| family", fam, "n", n,
    "quadrature", quad, "\n")
set.seed(20261021)
ng <- n / 6
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng), each = 6)))
u <- rnorm(ng, 0, 0.6)[d$g]
if (fam == "Beta") {
  m <- plogis(0.2 + 0.4 * d$x + u); d$y <- rbeta(n, m * 8, (1 - m) * 8)
  d$cc <- ifelse(d$y > 0.8, "right", "none"); d$y <- pmin(d$y, 0.8)
  F <- Beta()
} else if (fam == "gaussian") {
  d$y <- 0.2 + 0.4 * d$x + u + rnorm(n)
  d$cc <- ifelse(d$y > 1, "right", "none"); d$y <- pmin(d$y, 1)
  F <- gaussian()
} else if (fam == "poisson") {
  d$y <- rpois(n, exp(1 + 0.3 * d$x + u))
  d$cc <- ifelse(d$y >= 5, "right", "none"); d$y <- pmin(d$y, 5)
  F <- poisson()
}
t0 <- proc.time()[[3]]
f <- suppressWarnings(frm(bf(y | cens(cc) ~ x + (1 | g)), family = F,
                          data = d, quadrature = quad))
cat(sprintf("logLik %.6f | %.1f s\n", as.numeric(logLik(f)),
            proc.time()[[3]] - t0))
cat("END\n")
