# Fix round 1: the tape cost of the CDF pairs, per row: the time of one
# value-and-gradient sweep over 20,000 rows, against RTMB::ppois(), so
# that a family whose CDF a fit calls on every row (trunc(), or a draw
# of importance sampling) is priced.
#   Rscript dev/cdfs-fix1-cost.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), "\n")
n <- 20000
set.seed(1)
q <- stats::rpois(n, 3)
x0 <- log(3)
cost <- function(lab, f) {
  tp <- MakeTape(function(p) sum(f(p)), x0)
  t0 <- proc.time()[["elapsed"]]
  for (i in 1:20) tp$jacobian(x0)
  cat(sprintf("%-40s %8.2f ms per sweep\n", lab,
              1000 * (proc.time()[["elapsed"]] - t0) / 20))
}
cost("RTMB::ppois", function(p) log(RTMB::ppois(q, exp(p) + 0 * q)))
cost("poisson_cdf_log lS", function(p) ns$poisson_cdf_log(q, exp(p))$lS)
cost("poisson_cdf_log lF", function(p) ns$poisson_cdf_log(q, exp(p))$lF)
cost("log_igamma_two at a = q + 1", function(p) {
  ns$log_igamma_two(exp(p) + 0 * q, q + 1)$lP
})
cost("pnorm both tails", function(p) {
  z <- (q - exp(p)) / 2
  RTMB::pnorm(z, log.p = TRUE) + RTMB::pnorm(z, lower.tail = FALSE, log.p = TRUE)
})
cost("nbinom_cdf_log lS", function(p) ns$nbinom_cdf_log(q, p + 0 * q, 0.5)$lS)
cost("gamma_cdf_log lS", function(p) ns$gamma_cdf_log(q + 0.5, exp(p), 3)$lS)
cat("END\n")
