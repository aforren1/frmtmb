# Lane famend, 5722 and 3580: the tape cost of a CDF pair per row, and
# whether RTMB::pbeta() would be a cheaper binomial or negative binomial
# pair for the one-step residuals: its cost per sweep, its accuracy
# against stats::pbinom(log.p = TRUE), and the OSA time of a binomial
# (1 | g) fit with the pair swapped in this process only.
#   Rscript dev/famend-cost.R <lib or "base">
a <- commandArgs(trailingOnly = TRUE)
base <- "C:/Users/adf44/source/r/rellib-r21"
ulib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (a[1] != "base") a[1], base, ulib))
suppressPackageStartupMessages({library(RTMB); library(frmtmb)})
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
n <- 20000
set.seed(1)
q <- stats::rpois(n, 3)
qb <- pmin(stats::rbinom(n, 20, 0.3), 19)
x0 <- log(3)
cost <- function(lab, f, rounds = 5L) {
  tp <- MakeTape(function(p) sum(f(p)), x0)
  best <- Inf
  for (r in seq_len(rounds)) {
    t0 <- proc.time()[["elapsed"]]
    for (i in 1:20) tp$jacobian(x0)
    best <- min(best, (proc.time()[["elapsed"]] - t0) / 20)
  }
  cat(sprintf("COST %-44s %8.2f ms per sweep of %d rows\n", lab,
              1000 * best, n))
  1000 * best
}
c0 <- cost("RTMB::ppois (control)", function(p) {
  log(RTMB::ppois(q, exp(p) + 0 * q))
})
invisible(cost("poisson_cdf_log both", function(p) {
  l <- ns$poisson_cdf_log(q, exp(p)); l$lF + l$lS
}))
invisible(cost("binomial_cdf_log both", function(p) {
  pq <- list(p = 1 / (1 + exp(-p)) + 0 * qb, q = 1 / (1 + exp(p)) + 0 * qb)
  l <- ns$binomial_cdf_log(qb, 20, pq); l$lF + l$lS
}))
invisible(cost("RTMB::pbeta both, logged", function(p) {
  pr <- 1 / (1 + exp(-p)) + 0 * qb
  log(RTMB::pbeta(1 - pr, 20 - qb, qb + 1)) +
    log(RTMB::pbeta(pr, qb + 1, 20 - qb))
}))
invisible(cost("nbinom_cdf_log both", function(p) {
  l <- ns$nbinom_cdf_log(q, p + 0 * q, 0.5); l$lF + l$lS
}))

# accuracy of the logged RTMB::pbeta pair on doubles, binomial(size, 0.5)
for (sz in c(20, 200, 2000)) {
  qq <- c(0, 1, 2, floor(sz / 10), floor(sz / 2), sz - 3, sz - 1)
  pr <- 0.5
  lf <- log(RTMB::pbeta(1 - pr, sz - qq, qq + 1))
  ref <- stats::pbinom(qq, sz, pr, log.p = TRUE)
  cat(sprintf("ACC pbeta size %d: max |log F - ref| %.3g, non-finite %d of %d\n",
              sz, max(abs(lf - ref)[is.finite(lf)]), sum(!is.finite(lf)),
              length(qq)))
}

# the OSA time of a binomial (1 | g) fit, the pair as built and swapped
set.seed(7)
G <- 40
nn <- 800
d <- data.frame(x = stats::rnorm(nn), g = factor(rep(seq_len(G), each = 20)))
d$y <- stats::rbinom(nn, 20, stats::plogis(-0.5 + 0.4 * d$x +
                                             stats::rnorm(G, 0, 0.5)[d$g]))
d$k <- 20L
fit <- frm(y | trials(k) ~ x + (1 | g), family = binomial(), data = d)
t0 <- proc.time()[["elapsed"]]
r1 <- residuals(fit, type = "osa")[, "Estimate"]
t1 <- proc.time()[["elapsed"]]
orig <- ns$binomial_cdf_log
swap <- function(q, size, pq) {
  n <- max(length(q), length(size), length(pq[["p"]]))
  q <- rep_len(q, n)
  size <- rep_len(size, n)
  below <- q < 0
  above <- q >= size
  inside <- !below & !above
  qc <- ifelse(inside, q, 0)
  sz <- ifelse(inside, size, 1)
  lp2 <- list(lF = log(RTMB::pbeta(ns$cdf_rows(pq[["q"]], n), sz - qc,
                                   qc + 1)),
              lS = log(RTMB::pbeta(ns$cdf_rows(pq[["p"]], n), qc + 1,
                                   sz - qc)))
  ns$cdf_edges(lp2, below, above)
}
utils::assignInNamespace("binomial_cdf_log", swap, "frmtmb")
t2 <- proc.time()[["elapsed"]]
r2 <- tryCatch(residuals(fit, type = "osa")[, "Estimate"],
               error = function(e) {
                 cat("swapped ERROR", conditionMessage(e), "\n")
                 rep(NA, nn)
               })
t3 <- proc.time()[["elapsed"]]
utils::assignInNamespace("binomial_cdf_log", orig, "frmtmb")
cat(sprintf(paste0("OSA binomial (1 | g), %d rows: as built %.1f s, ",
                   "pbeta pair %.1f s; max |r1 - r2| %.3g\n"), nn,
            t1 - t0, t3 - t2, max(abs(r1 - r2))))
cat("END\n")
