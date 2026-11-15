# Lane famend, 3580: would a Horner sum of the poisson upper tail be a
# cheaper tape than the 40 logspace_add() terms of poisson_cdf_log()?
# Tape size (operations) and the time of one value-and-gradient sweep,
# rows right-censored at 5, against RTMB::ppois() as the control.
#   Rscript dev/famend-cost-pois.R <lib or "base">
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
n <- 20000
q <- rep(5, n)
x0 <- log(3)
cost <- function(lab, f, rounds = 5L) {
  tp <- MakeTape(function(p) sum(f(p)), x0)
  best <- Inf
  for (r in seq_len(rounds)) {
    t0 <- proc.time()[["elapsed"]]
    for (i in 1:20) tp$jacobian(x0)
    best <- min(best, (proc.time()[["elapsed"]] - t0) / 20)
  }
  ops <- tryCatch(length(tp$print(silent = TRUE)), error = function(e) NA)
  cat(sprintf("COST %-36s %8.2f ms per sweep, value %.10g\n", lab,
              1000 * best, tp(x0)))
}
horner_lS <- function(q, m, K = 40L) {
  lm <- log(m)
  mc <- ad_min(m, q + K + 1)
  rr <- ad_min(mc / (q + K + 1), exp(-0.01))
  s <- 1 / (1 - rr)
  for (k in seq.int(K, 2L)) s <- 1 + mc / (q + k) * s
  (q + 1) * lm - m - lgamma(q + 2) + log(s)
}
ad_min <- ns$ad_min
cost("RTMB::ppois (control)", function(p) log(RTMB::ppois(q, exp(p) + 0 * q)))
cost("poisson_cdf_log lS", function(p) ns$poisson_cdf_log(q, exp(p) + 0 * q)$lS)
cost("Horner upper sum lS", function(p) horner_lS(q, exp(p) + 0 * q))
for (m in c(0.5, 3, 6, 20)) {
  cat(sprintf("ACC m %4.1f: horner %.15g ppois %.15g\n", m,
              horner_lS(5, m), ppois(5, m, lower.tail = FALSE, log.p = TRUE)))
}
