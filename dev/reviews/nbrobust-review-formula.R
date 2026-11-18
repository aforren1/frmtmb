# Review of lane nbrobust, check 1: nb_lpmf() against an independent
# high-precision reference. The reference value is the textbook NB2 mass
#   lgamma(y + n) - lgamma(n) - lgamma(y + 1) + n log(n / (n + mu))
#     + y log(mu / (n + mu))
# in Rmpfr at 400 bits, and its derivatives in log mu and log n are
# CENTRAL DIFFERENCES of that mpfr value (step 2^-100), so no analytic
# derivative formula is shared with the lane.
#   Rscript dev/reviews/nbrobust-review-formula.R
.libPaths(c("C:/Users/adf44/source/r/wt-nbrobust-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(Rmpfr)
})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
nb_lpmf <- frmtmb:::nb_lpmf
lsd <- frmtmb:::lgamma_shift_diff
eps <- .Machine$double.eps
P <- 400
ref <- function(y, lmu, lsh) {
  f <- function(a, b) {
    mu <- exp(a)
    n <- exp(b)
    lgamma(y + n) - lgamma(n) - lgamma(mpfr(y + 1, P)) +
      n * (b - log(n + mu)) + y * (a - log(n + mu))
  }
  a <- mpfr(lmu, P)
  b <- mpfr(lsh, P)
  h <- mpfr(2, P)^-100
  v <- f(a, b)
  da <- (f(a + h, b) - f(a - h, b)) / (2 * h)
  db <- (f(a, b + h) - f(a, b - h)) / (2 * h)
  c(v = as.numeric(v), dmu = as.numeric(da), dsh = as.numeric(db))
}
# arms on the tape; "vec" gives a shape per row (nb_rising's loop branch),
# "osa" the lgamma_shift_diff() form the osa route takes for every row,
# "robust" the 0.83.0 call
arm_fun <- function(arm, y) {
  switch(arm,
    scalar = function(p) nb_lpmf(c(y, y), p[1], p[2])[1],
    vec = function(p) nb_lpmf(c(y, y), c(p[1], p[1]), c(p[2], p[2]))[1],
    osa = function(p) {
      lsh <- p[2]
      lmu <- p[1]
      lsd(exp(lsh), y) - y * lsh + y * lmu -
        (exp(lsh) + y) * RTMB::logspace_add(0 * lmu, lmu - lsh) -
        lgamma(y + 1)
    },
    robust = function(p) RTMB::dnbinom_robust(y, p[1], 2 * p[1] - p[2],
                                              log = TRUE))
}
ys <- c(0, 1, 2, 3, 7, 23, 24, 25, 26, 60, 500, 1e4, 1e5)
lmus <- c(-30, -5, -1, 0, 1.5, 5, 9.2, 11.5)
lshs <- c(-10, -3, 0, 2, 5, 10, 15, 18, 20, 22, 24, 26, 28, 30, 35, 40, 50)
rows <- list()
for (y in ys) for (lmu in lmus) for (lsh in lshs) {
  r <- ref(y, lmu, lsh)
  sd <- stats::dnbinom(y, size = exp(lsh), mu = exp(lmu), log = TRUE)
  for (arm in c("scalar", "vec", "osa", "robust")) {
    f <- arm_fun(arm, y)
    tp <- RTMB::MakeTape(f, c(lmu, lsh))
    v <- tp(c(lmu, lsh))
    g <- as.numeric(tp$jacobian(c(lmu, lsh)))
    # re-evaluate the tape at another point and back, so a value fixed
    # at taping would show
    v2 <- tp(c(lmu + 0.5, lsh + 1))
    r2 <- ref(y, lmu + 0.5, lsh + 1)[["v"]]
    unit <- eps * (y + exp(lmu) + 1)
    rows[[length(rows) + 1L]] <- data.frame(
      arm = arm, y = y, lmu = lmu, lsh = lsh, v = v, vref = r[["v"]],
      vstats = sd,
      relv = abs(v - r[["v"]]) / max(abs(r[["v"]]), 1e-300),
      absv = abs(v - r[["v"]]),
      absv_retape = abs(v2 - r2),
      dmu = g[1], dmu_ref = r[["dmu"]], dsh = g[2], dsh_ref = r[["dsh"]],
      udmu = abs(g[1] - r[["dmu"]]) / unit,
      udsh = abs(g[2] - r[["dsh"]]) / unit)
  }
}
tab <- do.call(rbind, rows)
write.table(tab, "dev/reviews/nbrobust-review-formula.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
band <- cut(tab$lsh, c(-Inf, 10, 18, 24, 30, Inf),
            labels = c("<=10", "10-18", "18-24", "24-30", ">30"))
cat("\nmax over the grid, by arm and log-shape band\n")
for (a in unique(tab$arm)) for (b in levels(band)) {
  s <- tab[tab$arm == a & band == b, ]
  cat(sprintf(paste("%-7s lsh %-6s n %4d  max units dsh %10.3g  dmu %8.3g",
                    " max abs val %9.3g  max rel val %9.3g  retape %9.3g\n"),
              a, b, nrow(s), max(s$udsh), max(s$udmu), max(s$absv),
              max(s$relv), max(s$absv_retape)))
}
s <- tab[tab$arm == "scalar", ]
cat("\nscalar arm, worst dsh units:\n")
print(head(s[order(-s$udsh), c("y", "lmu", "lsh", "dsh", "dsh_ref", "udsh")],
           5))
cat("scalar arm, worst abs value error:\n")
print(head(s[order(-s$absv), c("y", "lmu", "lsh", "v", "vref", "absv")], 5))
cat("\nstats::dnbinom against the mpfr reference (abs), max:",
    format(max(abs(s$vstats - s$vref)), digits = 3), "\n")
cat("nb_lpmf against stats::dnbinom, max abs", format(max(abs(s$v - s$vstats)),
    digits = 3), " max rel", format(max(abs(s$v - s$vstats) /
    pmax(abs(s$vstats), 1e-300)), digits = 3), "\n")
# y = 24 on both branches, and y = 25, value and derivative
cat("\nswitch at 24: the sum form against lgamma_shift_diff() at the same y\n")
for (y in c(2, 23, 24, 25)) for (lsh in c(-3, 0, 5, 20, 35)) {
  f1 <- function(p) {
    j <- seq_len(max(y - 1, 0))
    if (length(j)) sum(log1p(j * exp(-p))) else 0 * p
  }
  f2 <- function(p) lsd(exp(p), y) - y * p
  t1 <- RTMB::MakeTape(f1, lsh)
  t2 <- RTMB::MakeTape(f2, lsh)
  cat(sprintf("y %2d lsh %3g  value diff %9.3g  deriv diff %9.3g (unit %9.3g)\n",
              y, lsh, t1(lsh) - t2(lsh),
              as.numeric(t1$jacobian(lsh)) - as.numeric(t2$jacobian(lsh)),
              eps * (y + 1)))
}
# what nb_rising() does with a count the validator would refuse, and NA
cat("\nnon-integer y = 2.5 at lsh 0, lmu 0: ", nb_lpmf(2.5, 0, 0),
    " stats:", suppressWarnings(stats::dnbinom(2.5, size = 1, mu = 1,
                                               log = TRUE)),
    " gamma-extension ref:", ref(2.5, 0, 0)[["v"]], "\n")
cat("NA y: ")
print(tryCatch(nb_lpmf(c(1, NA, 3), 0, 0), error = function(e) conditionMessage(e)))
cat("NA y, robust: ")
print(RTMB::dnbinom_robust(c(1, NA, 3), 0, 0, log = TRUE))
