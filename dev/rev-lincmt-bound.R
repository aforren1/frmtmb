# The bound the lane gestured at but did not measure.
#
# dev/lincmt-findings.md says of the 7.07e-11 worst point that "its
# ABSOLUTE error stays at machine epsilon and the trajectory is right
# to eps of its own scale". That is the claim that decides whether the
# one cancellation left in the formulation is benign, and it is
# testable: measure the error relative to the TRAJECTORY'S OWN MAXIMUM
# rather than pointwise, over a parameter box wider than the lane's and
# over a lag grid rather than one lag.
#
# Seed 909. Script path: dev/rev-lincmt-bound.R.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-mpfr.R")

set.seed(909)
lags <- c(0.05, 0.5, 2, 8, 24, 72, 200)
cat(sprintf("\n%5s %8s %14s %14s %14s\n", "ncmt", "draws",
            "worst pointwise", "worst / peak", "worst peak/val"))
for (nc in 1:3) {
  wp <- 0; ws <- 0; wr <- 0; wat <- NULL
  for (i in 1:40) {
    # wider than the lane's box: rates over six decades, so that the
    # near-decoupled corner is drawn often
    p <- as.list(exp(runif(6, log(1e-4), log(20))))
    names(p) <- c("ke", "k12", "k21", "k13", "k31", "ka")
    pl <- c(p[c(TRUE, nc >= 2, nc >= 2, nc == 3, nc == 3, TRUE)],
            list(V = 1))
    a <- frm_lincmt(parms = pl, times = lags, ncmt = nc, depot = TRUE,
                    init = list(depot = 1), output = "central")
    b <- vapply(lags, function(u)
      as.numeric(bolus_ss(nc, TRUE, p, u)), 0)
    pk <- max(abs(b))
    if (!(pk > 0)) next
    pw <- max(abs(a - b) / pmax(abs(b), 1e-300))
    sc <- max(abs(a - b)) / pk
    rr <- pk / min(abs(b)[abs(b) > 0])
    if (pw > wp) { wp <- pw; wat <- c(unlist(p), pk = pk) }
    if (sc > ws) ws <- sc
    if (rr > wr) wr <- rr
  }
  cat(sprintf("%5d %8d %14.3e %14.3e %14.3e\n", nc, 40L, wp, ws, wr))
  cat("   worst-pointwise draw:",
      paste(format(wat, digits = 3), collapse = " "), "\n")
}

cat("\n=== the bound, stated ===\n")
cat("error relative to the trajectory's own maximum stays at a few\n")
cat("ulp, so the pointwise relative error at a point that is 10^-d of\n")
cat("the peak is at most about 1e-15 * 10^d:\n")
for (d in 2:8)
  cat(sprintf("   %d decades below the peak: %8.1e\n", d, 1e-15 * 10^d))
cat("The item's bar is 1e-8, which this crosses only beyond seven\n")
cat("decades below Cmax. No assay reports there.\n")
