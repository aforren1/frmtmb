# Attack 1 and 2: how far does the one cancellation the lane left in
# the two-compartment disposition actually go, and is the 7.07e-11 the
# tail of a bounded distribution or the visible edge of a region?
#
# Run: Rscript rev-lincmt-cancel.R
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-mpfr.R")

lin <- function(p, u, ncmt, depot = TRUE) {
  frm_lincmt(parms = c(p, list(V = 1)), times = u, ncmt = ncmt,
             depot = depot,
             init = if (depot) list(depot = 1) else list(central = 1),
             output = "central")
}

cat("\n=== A. the two references agree with each other ===\n")
set.seed(4041)
wa <- 0
for (i in 1:12) {
  p <- as.list(exp(runif(6, log(1e-3), log(5))))
  names(p) <- c("ke", "k12", "k21", "k13", "k31", "ka")
  u <- exp(runif(1, log(0.05), log(72)))
  for (nc in 1:3) {
    a <- bolus_ss(nc, TRUE, p, u)
    b <- bolus_pf(nc, TRUE, p, u)
    r <- as.numeric(abs(a - b) / abs(b))
    if (is.finite(r) && r > wa) wa <- r
  }
}
cat("worst scaling-and-squaring vs partial-fraction, 300 bits:",
    format(wa), "\n")

cat("\n=== B. the lane's random sweep, seed 11, nc = 2, reproduced ===\n")
set.seed(11)
for (nc in 1:3) {
  worst <- 0; wat <- NULL
  for (i in 1:60) {
    p <- as.list(exp(runif(6, log(1e-3), log(5))))
    names(p) <- c("ke", "k12", "k21", "k13", "k31", "ka")
    u <- exp(runif(1, log(0.05), log(72)))
    pp <- p[c(TRUE, nc >= 2, nc >= 2, nc == 3, nc == 3, TRUE)]
    a <- lin(pp, u, nc)
    b <- bolus_ss(nc, TRUE, p, u)
    r <- relerr(a, b)
    if (is.finite(r) && r > worst) { worst <- r; wat <- c(unlist(p), u) }
  }
  cat(nc, "cmt + depot: worst rel", format(worst, digits = 4), "at",
      paste(format(wat, digits = 4), collapse = " "), "\n")
}

cat("\n=== C. the named worst point, decomposed ===\n")
p2 <- list(ke = 2.393, k12 = 0.001689, k21 = 0.004481, ka = 0.4100)
u2 <- 34.88
g <- (p2$ke + p2$k12 - p2$k21) / 2
dd <- sqrt(g * g + p2$k12 * p2$k21)
c1 <- (dd + g) / (2 * dd)
lam1 <- (p2$ke + p2$k12 + p2$k21) / 2 + dd
lam2 <- p2$ke * p2$k21 / lam1
cat(sprintf("g %g  dd %g  c_small %g  k12k21/g^2 %g\n",
            g, dd, 1 - c1, p2$k12 * p2$k21 / (g * g)))
cat(sprintf("lam1 %g  lam2 %g  exp(-lam1 u)/exp(-lam2 u) %g\n",
            lam1, lam2, exp(-(lam1 - lam2) * u2)))
a <- lin(p2, u2, 2L); b <- bolus_ss(2L, TRUE, p2, u2)
cat(sprintf("pointwise rel err %.3e ; predicted eps/(4 c_small) %.3e\n",
            relerr(a, b), 2.22e-16 / (4 * (1 - c1))))

cat("\n=== D. push it: k12 k21 / g^2 driven down ===\n")
cat("2 cmt + depot, ke = 2, ka = 0.41, k12 = k21 = r.\n")
cat("`peak` is the maximum of the closed form over a 0.05 to 200 grid;",
    "\n`errS` is the absolute error at the reported lag divided by that",
    "peak.\n\n")
cat(sprintf("%10s %12s %8s %14s %12s %12s\n", "r", "k12k21/g^2",
            "u", "value", "relerr", "errS"))
grid <- exp(seq(log(0.05), log(200), length.out = 400))
for (r in 10^-(1:8)) {
  p <- list(ke = 2, k12 = r, k21 = r, ka = 0.41)
  gg <- (p$ke + p$k12 - p$k21) / 2
  tr <- lin(p, grid, 2L)
  peak <- max(abs(tr))
  for (u in c(10, 40, 120)) {
    a <- lin(p, u, 2L)
    b <- bolus_ss(2L, TRUE, p, u)
    cat(sprintf("%10.1e %12.2e %8g %14.6e %12.3e %12.3e\n",
                r, p$k12 * p$k21 / (gg * gg), u, a, relerr(a, b),
                as.numeric(abs(mp(a) - b)) / peak))
  }
}

cat("\n=== E. the same sweep on the TEXTBOOK spelling, for contrast ==\n")
cat("textbook 2 cmt + depot central amount, double precision.\n")
tb <- function(p, u) {
  b <- p$ke + p$k12 + p$k21
  d <- sqrt(b * b - 4 * p$ke * p$k21)
  l1 <- (b + d) / 2; l2 <- (b - d) / 2
  c1 <- (l1 - p$k21) / (l1 - l2); c2 <- (l2 - p$k21) / (l2 - l1)
  ka <- p$ka
  ka * (c1 * (exp(-l1 * u) - exp(-ka * u)) / (ka - l1) +
        c2 * (exp(-l2 * u) - exp(-ka * u)) / (ka - l2))
}
cat(sprintf("%10s %8s %12s %12s\n", "r", "u", "shipped", "textbook"))
for (r in 10^-(1:8)) {
  p <- list(ke = 2, k12 = r, k21 = r, ka = 0.41)
  for (u in c(40, 120)) {
    b <- bolus_ss(2L, TRUE, p, u)
    cat(sprintf("%10.1e %8g %12.3e %12.3e\n", r, u,
                relerr(lin(p, u, 2L), b), relerr(tb(p, u), b)))
  }
}

cat("\n=== F. near ka == ke, shipped vs textbook, 1 cmt ===\n")
cat(sprintf("%10s %12s %12s\n", "ka - ke", "shipped", "textbook"))
tb1 <- function(p, u)
  p$ka * (exp(-p$ke * u) - exp(-p$ka * u)) / (p$ka - p$ke)
for (d in c(0, 10^-(16:2))) {
  p <- list(ke = 0.2, ka = 0.2 + d)
  u <- 6
  b <- bolus_ss(1L, TRUE, p, u)
  cat(sprintf("%10.1e %12.3e %12.3e\n", d, relerr(lin(p, u, 1L), b),
              relerr(tb1(p, u), b)))
}
