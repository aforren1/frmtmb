# Round 2, item 1, part six, and the number that decides the severity.
#
# dev/rev-lincmt-f6e.R found that on a ten-decade rate box the VALUE is
# wrong by up to 1.3e-04 of the trajectory's own maximum AT DRAWS WHERE
# THE GRADIENT IS FINITE, that is, silently. The NaN is a symptom of a
# neighbourhood, not a knife edge.
#
# The severity turns on where that neighbourhood starts. This measures
# the worst value error against the 300-bit reference as a function of
# the one thing a user can see, the spread of the rate constants, for
# two and for three compartments. Two compartments solve a quadratic
# and never call acos(), so they are the control.
#
# Script path: dev/rev-lincmt-f6f.R. Seed 1234.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-mpfr.R")

lags <- c(0.2, 2, 12, 96)
err_of <- function(nc, p) {
  pl <- c(p[c("ke", "k12", "k21", "k13", "k31",
              "ka")[c(TRUE, nc >= 2, nc >= 2, nc == 3, nc == 3,
                      TRUE)]], list(V = 1))
  a <- tryCatch(frm_lincmt(parms = pl, times = lags, ncmt = nc,
                           depot = TRUE, init = list(depot = 1),
                           output = "central"),
                error = function(e) NULL)
  if (is.null(a) || anyNA(a)) return(NA_real_)
  b <- vapply(lags, function(u) as.numeric(bolus_ss(nc, TRUE, p, u)),
              0)
  pk <- max(abs(b))
  if (!(pk > 0) || anyNA(b)) return(NA_real_)
  max(abs(a - b)) / pk
}

set.seed(1234)
N <- 40L
cat("\n=== worst value error / peak, against the rate spread ===\n")
cat(N, "draws per cell, rates log-uniform over a window of the given",
    "\nwidth centred on 0.2 per hour. `worst` is the largest error",
    "\nrelative to the trajectory's own maximum.\n\n")
cat(sprintf("%8s %14s %14s %14s %14s\n", "decades", "2 cmt worst",
            "2 cmt median", "3 cmt worst", "3 cmt median"))
for (dec in c(2, 4, 6, 8, 10)) {
  out <- list()
  for (nc in 2:3) {
    e <- numeric(0)
    for (i in seq_len(N)) {
      p <- as.list(0.2 * 10^runif(6, -dec / 2, dec / 2))
      names(p) <- c("ke", "k12", "k21", "k13", "k31", "ka")
      v <- err_of(nc, p)
      if (is.finite(v)) e <- c(e, v)
    }
    out[[nc - 1L]] <- e
  }
  cat(sprintf("%8d %14.3e %14.3e %14.3e %14.3e\n", dec,
              max(out[[1L]]), median(out[[1L]]),
              max(out[[2L]]), median(out[[2L]])))
}

cat("\n=== the same, restricted to a PK-shaped model ===\n")
cat("A real three-compartment PK model is not six independent draws.",
    "\nka and ke sit within a decade or two of each other, and the",
    "\nperipheral rates are the ones that can be small. This draws ka",
    "\nand ke over two decades around 0.2 and lets the four peripheral",
    "\nrates span the given width downward, which is what an",
    "\noptimizer does to a compartment the data do not support.\n\n")
set.seed(1234)
cat(sprintf("%18s %14s %14s %10s\n", "peripheral floor", "worst",
            "median", "NaN grad"))
gradna <- function(p) {
  f <- function(th) sum(frm_lincmt(
    parms = list(ke = exp(th[1]), k12 = exp(th[2]), k21 = exp(th[3]),
                 k13 = exp(th[4]), k31 = exp(th[5]), ka = exp(th[6]),
                 V = 1),
    times = lags, ncmt = 3, depot = TRUE, init = list(depot = 1),
    output = "central"))
  x <- log(unlist(p[c("ke", "k12", "k21", "k13", "k31", "ka")]))
  !all(is.finite(as.numeric(MakeTape(f, x)$jacobian(x))))
}
for (flo in c(-4, -8, -10, -12)) {
  e <- numeric(0); nn <- 0L
  for (i in seq_len(N)) {
    p <- list(ke = 0.2 * 10^runif(1, -1, 1),
              k12 = 10^runif(1, flo, 0.5), k21 = 10^runif(1, flo, 0.5),
              k13 = 10^runif(1, flo, 0.5), k31 = 10^runif(1, flo, 0.5),
              ka = 0.2 * 10^runif(1, -1, 1))
    v <- err_of(3L, p)
    if (is.finite(v)) e <- c(e, v) else nn <- nn + 1L
    if (gradna(p)) nn <- nn + 0L
  }
  cat(sprintf("%18s %14.3e %14.3e %10d\n", paste0("1e", flo),
              max(e), median(e), nn))
}
