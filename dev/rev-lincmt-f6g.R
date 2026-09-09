# The domain question, answered with a cheap reference.
#
# dev/rev-lincmt-f6f.R asks where the three-compartment value starts to
# degrade with the spread of the rate constants, against the 300-bit
# reference. At ten decades of spread the scaling and squaring needs so
# many squarings that the sweep does not finish in a reasonable time.
#
# It does not need 300 bits. frm_ode() at atol = rtol = 1e-12 is an
# independent implementation whose own error is about 1e-12 of the
# trajectory's scale, which is four decades below the 1e-08 the item
# asks about and eight below the 1.3e-04 dev/rev-lincmt-f6e.R found.
# So it can settle where the degradation begins even though it cannot
# settle the last few digits.
#
# Script path: dev/rev-lincmt-f6g.R. Seed 1234, the same draws
# dev/rev-lincmt-f6f.R uses.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")

lags <- c(0.2, 2, 12, 96)
d3 <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[[6]] * y[1],
         p[[6]] * y[1] - (p[[1]] + p[[2]] + p[[4]]) * y[2] +
           p[[3]] * y[3] + p[[5]] * y[4],
         p[[2]] * y[2] - p[[3]] * y[3],
         p[[4]] * y[2] - p[[5]] * y[4]))
}
d2 <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[[6]] * y[1],
         p[[6]] * y[1] - (p[[1]] + p[[2]]) * y[2] + p[[3]] * y[3],
         p[[2]] * y[2] - p[[3]] * y[3]))
}
err_of <- function(nc, p) {
  pl <- c(p[c("ke", "k12", "k21", "k13", "k31",
              "ka")[c(TRUE, nc >= 2, nc >= 2, nc == 3, nc == 3,
                      TRUE)]], list(V = 1))
  a <- tryCatch(frm_lincmt(parms = pl, times = lags, ncmt = nc,
                           depot = TRUE, init = list(depot = 1),
                           output = "central"),
                error = function(e) NULL)
  if (is.null(a) || anyNA(a)) return(NA_real_)
  pv <- list(p$ke, p$k12, p$k21, p$k13, p$k31, p$ka)
  st <- if (nc == 3L) c("depot", "central", "p1", "p2") else
    c("depot", "central", "p1")
  # the same unit bolus into the depot the closed form is given
  i0 <- as.list(rep(0, nc + 1L)); i0[[1L]] <- 1
  one <- function(tol) tryCatch(
    frm_ode(if (nc == 3L) d3 else d2,
            init = i0, times = lags,
            parms = pv, states = st, output = "central",
            atol = tol, rtol = tol), error = function(e) NULL)
  b <- one(1e-12)
  if (is.null(b) || anyNA(b)) return(NA_real_)
  l <- one(1e-9)
  if (is.null(l) || anyNA(l)) return(NA_real_)
  pk <- max(abs(b))
  if (!(pk > 0)) return(NA_real_)
  # report the disagreement AND what the solver's own tolerance costs,
  # so that a number can be read as the closed form's only when it is
  # well above the solver's own spread
  c(max(abs(a - b)) / pk, max(abs(l - b)) / pk)
}

set.seed(1234)
N <- 60L
cat("\n=== value error against the rate spread ===\n")
cat(N, "draws per cell, rates log-uniform over a window of the given",
    "\nwidth centred on 0.2 per hour, one bolus into the depot, four",
    "\nlags. `lincmt` is the closed form against frm_ode() at 1e-12;",
    "\n`solver` is what frm_ode()'s own tolerance costs on the same",
    "\ndraws, which is the floor below which the first column means",
    "\nnothing.\n\n")
cat(sprintf("%8s %6s %12s %12s %12s %12s\n", "decades", "ncmt",
            "lincmt max", "lincmt med", "solver max", "n"))
for (dec in c(2, 4, 6, 8, 10, 12)) {
  for (nc in 2:3) {
    e <- numeric(0); sv <- numeric(0)
    for (i in seq_len(N)) {
      p <- as.list(0.2 * 10^runif(6, -dec / 2, dec / 2))
      names(p) <- c("ke", "k12", "k21", "k13", "k31", "ka")
      v <- suppressWarnings(err_of(nc, p))
      if (length(v) == 2L && all(is.finite(v))) {
        e <- c(e, v[[1L]]); sv <- c(sv, v[[2L]])
      }
    }
    if (!length(e)) { cat(sprintf("%8d %6d %12s\n", dec, nc, "-")); next }
    cat(sprintf("%8d %6d %12.3e %12.3e %12.3e %12d\n", dec, nc,
                max(e), median(e), max(sv), length(e)))
  }
}

cat("\n=== a PK-shaped model, the peripheral rates driven down ===\n")
cat("ka and ke over two decades around 0.2, the four peripheral rates",
    "\nspanning from the given floor up to about 3, which is what an",
    "\noptimizer does to a compartment the data do not support.\n\n")
set.seed(1234)
cat(sprintf("%12s %12s %12s %12s %8s\n", "floor", "lincmt max",
            "lincmt med", "solver max", "n"))
for (flo in c(-2, -4, -6, -8, -10, -12)) {
  e <- numeric(0); sv <- numeric(0)
  for (i in seq_len(N)) {
    p <- list(ke = 0.2 * 10^runif(1, -1, 1),
              k12 = 10^runif(1, flo, 0.5), k21 = 10^runif(1, flo, 0.5),
              k13 = 10^runif(1, flo, 0.5), k31 = 10^runif(1, flo, 0.5),
              ka = 0.2 * 10^runif(1, -1, 1))
    v <- suppressWarnings(err_of(3L, p))
    if (length(v) == 2L && all(is.finite(v))) {
      e <- c(e, v[[1L]]); sv <- c(sv, v[[2L]])
    }
  }
  if (!length(e)) { cat(sprintf("%12s %12s\n", paste0("1e", flo), "-"))
    next }
  cat(sprintf("%12s %12.3e %12.3e %12.3e %8d\n", paste0("1e", flo),
              max(e), median(e), max(sv), length(e)))
}
