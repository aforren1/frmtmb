# Attack: an eigenvalue coalescence the lane's sweeps do NOT cover.
#
# lincmt_disp()'s three-compartment coefficient is
#   cf(x, y, z) = (x - k21)(x - k31) / ((x - y)(x - z) + 1e-150),
# which divides by a DIFFERENCE OF EIGENVALUES. For a mammillary model
# with k12 > 0 and k13 > 0 the poles and the zeros strictly interlace,
# so a coalescence forces the numerator to vanish too and the quotient
# stays finite. That protection is lost when k13 (or k12) is zero or
# nearly zero: the third compartment decouples, k31 becomes an
# eigenvalue whose numerator factor cancels, and k31 is free to collide
# with a root of the remaining quadratic.
#
# This script asks what frm_lincmt() returns there.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-mpfr.R")

lin <- function(p, u, ncmt, depot = TRUE)
  frm_lincmt(parms = c(p, list(V = 1)), times = u, ncmt = ncmt,
             depot = depot,
             init = if (depot) list(depot = 1) else list(central = 1),
             output = "central")

# the two disposition eigenvalues that survive when k13 = 0
quad_root <- function(ke, k12, k21) {
  b <- ke + k12 + k21
  d <- sqrt(b * b - 4 * ke * k21)
  c(hi = (b + d) / 2, lo = (b - d) / 2)
}

KE <- 0.2; K12 <- 0.4; K21 <- 0.1
qr <- quad_root(KE, K12, K21)
cat(sprintf("\nreduced eigenvalues with k13 = 0: %.15g and %.15g\n",
            qr[["hi"]], qr[["lo"]]))

cat("\n=== A. k13 exactly 0, k31 exactly on the slow root ===\n")
p <- list(ke = KE, k12 = K12, k21 = K21, k13 = 0, k31 = qr[["lo"]],
          ka = 1.1)
grid <- c(0.5, 1, 2, 4, 8, 12, 24, 48)
a <- lin(p, grid, 3L)
b <- vapply(grid, function(u) as.numeric(bolus_ss(3L, TRUE, p, u)), 0)
cat(sprintf("%6s %16s %16s %12s\n", "u", "frm_lincmt", "300-bit",
            "rel"))
for (i in seq_along(grid))
  cat(sprintf("%6g %16.8e %16.8e %12.3e\n", grid[[i]], a[[i]], b[[i]],
              abs(a[[i]] - b[[i]]) / abs(b[[i]])))
cat("error at the peak, relative to the peak:",
    format(max(abs(a - b)) / max(abs(b)), digits = 4), "\n")

cat("\n=== B. does frm_ode() agree with the 300-bit reference? ===\n")
dyn <- function(t, y, pr) list(c(
  -pr[[6]] * y[1],
  pr[[6]] * y[1] - (pr[[1]] + pr[[2]] + pr[[4]]) * y[2] +
    pr[[3]] * y[3] + pr[[5]] * y[4],
  pr[[2]] * y[2] - pr[[3]] * y[3],
  pr[[4]] * y[2] - pr[[5]] * y[4]))
o <- frm_ode(dyn, init = list(1, 0, 0, 0), times = grid,
             parms = list(p$ke, p$k12, p$k21, p$k13, p$k31, p$ka),
             states = c("depot", "central", "p1", "p2"),
             output = "central", atol = 1e-12, rtol = 1e-12)
cat(sprintf("%6s %16s %16s %12s\n", "u", "frm_ode", "300-bit", "rel"))
for (i in seq_along(grid))
  cat(sprintf("%6g %16.8e %16.8e %12.3e\n", grid[[i]], o[[i]], b[[i]],
              abs(o[[i]] - b[[i]]) / abs(b[[i]])))

cat("\n=== C. k13 driven to zero, k31 pinned on the slow root ===\n")
cat(sprintf("%10s %10s %14s %14s %12s %12s\n", "k13", "split", "lincmt",
            "300-bit", "rel", "rel/peak"))
for (k13 in c(10^-(1:14), 0)) {
  p <- list(ke = KE, k12 = K12, k21 = K21, k13 = k13,
            k31 = qr[["lo"]], ka = 1.1)
  a <- lin(p, grid, 3L)
  b <- vapply(grid, function(u)
    as.numeric(bolus_ss(3L, TRUE, p, u)), 0)
  # how far apart the two colliding eigenvalues are, in double
  d3 <- frmtmb.ode:::lincmt_disp(3L, p$ke, p$k12, p$k21, p$k13, p$k31)
  lv <- sort(vapply(d3$lam, as.numeric, 0))
  spl <- min(diff(lv))
  j <- which.max(abs(a - b))
  cat(sprintf("%10.1e %10.2e %14.6e %14.6e %12.3e %12.3e\n", k13, spl,
              a[[j]], b[[j]], abs(a[[j]] - b[[j]]) / abs(b[[j]]),
              max(abs(a - b)) / max(abs(b))))
}

cat("\n=== D. k31 swept across the resonance, k13 = 1e-6 ===\n")
cat(sprintf("%14s %14s %12s %12s\n", "k31", "worst rel", "rel/peak",
            "split"))
for (fac in c(0.5, 0.9, 0.99, 0.999, 1, 1.001, 1.01, 1.1, 2)) {
  p <- list(ke = KE, k12 = K12, k21 = K21, k13 = 1e-6,
            k31 = qr[["lo"]] * fac, ka = 1.1)
  a <- lin(p, grid, 3L)
  b <- vapply(grid, function(u)
    as.numeric(bolus_ss(3L, TRUE, p, u)), 0)
  d3 <- frmtmb.ode:::lincmt_disp(3L, p$ke, p$k12, p$k21, p$k13, p$k31)
  lv <- sort(vapply(d3$lam, as.numeric, 0))
  cat(sprintf("%14.8f %14.3e %12.3e %12.3e\n", p$k31,
              max(abs(a - b) / abs(b)), max(abs(a - b)) / max(abs(b)),
              min(diff(lv))))
}

cat("\n=== E. the same collision reached from k12 instead ===\n")
qr2 <- quad_root(KE, 0.4, 0.1)
cat(sprintf("%10s %14s %12s\n", "k12", "worst rel", "rel/peak"))
for (k12 in c(10^-(1:12), 0)) {
  # k13/k31 carry the model; k12 -> 0 decouples peripheral 1 and lets
  # k21 collide with a root of the (ke, k13, k31) quadratic
  qq <- quad_root(KE, 0.4, 0.1)
  p <- list(ke = KE, k12 = k12, k21 = qq[["lo"]], k13 = 0.4,
            k31 = 0.1, ka = 1.1)
  a <- lin(p, grid, 3L)
  b <- vapply(grid, function(u)
    as.numeric(bolus_ss(3L, TRUE, p, u)), 0)
  cat(sprintf("%10.1e %14.3e %12.3e\n", k12,
              max(abs(a - b) / abs(b)), max(abs(a - b)) / max(abs(b))))
}

cat("\n=== F. rate constants that underflow to zero ===\n")
show <- function(tag, ncmt, p, depot = TRUE) {
  v <- tryCatch(lin(p, c(1, 6, 24), ncmt, depot),
                error = function(e) conditionMessage(e))
  cat(sprintf("%-34s %s\n", tag,
              if (is.character(v)) paste("ERROR:", substr(v, 1, 40))
              else paste(format(as.numeric(v), digits = 6),
                         collapse = "  ")))
}
show("1 cmt, ke = 0", 1L, list(ke = 0, ka = 1.1))
show("2 cmt, every rate 0", 2L,
     list(ke = 0, k12 = 0, k21 = 0, ka = 1.1))
show("3 cmt, every rate 0", 3L,
     list(ke = 0, k12 = 0, k21 = 0, k13 = 0, k31 = 0, ka = 1.1))
show("3 cmt, ke = 0 only", 3L,
     list(ke = 0, k12 = 0.4, k21 = 0.1, k13 = 0.4, k31 = 0.1,
          ka = 1.1))
show("3 cmt, k21 = k31 = 0", 3L,
     list(ke = 0.2, k12 = 0.4, k21 = 0, k13 = 0.4, k31 = 0,
          ka = 1.1))
show("3 cmt, k12 = k13 = 0", 3L,
     list(ke = 0.2, k12 = 0, k21 = 0.1, k13 = 0, k31 = 0.3,
          ka = 1.1))
show("2 cmt, ke = 0, k12 = k21 = 0.2", 2L,
     list(ke = 0, k12 = 0.2, k21 = 0.2, ka = 1.1))
