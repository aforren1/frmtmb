# F6 followed through: the review asked for the eigenvalue collision
# that the 72-case gradient sweep does not contain. Adding it found a
# boundary the review did not have.
#
# `lincmt_disp()` solves the cubic through `acos()`, whose derivative is
# infinite where a double root puts its argument at exactly one. This
# measures where that bites: how large the eigenvalue split has to be
# before the gradient is finite, and whether the VALUE is right on the
# other side of that line.
#
# Seed: none, every point is constructed. Run from this directory.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src2.R")

tt <- c(0.5, 1, 2, 4, 8, 12, 24, 48)
ev <- data.frame(time = 0, state = "depot", value = 100, ii = 8,
                 addl = 3L)
KE <- 0.2; K12 <- 0.4; K21 <- 0.1; K13 <- 1e-300
b <- KE + K12 + K21
# the slow root of the reduced quadratic; k13 at 1e-300 decouples the
# third compartment, so putting k31 here makes that root double
qlo <- (b - sqrt(b * b - 4 * KE * K21)) / 2
alt <- KE * K21 / ((b + sqrt(b * b - 4 * KE * K21)) / 2)
cat("the slow root, two algebraically identical spellings:\n  ",
    format(qlo, digits = 17), "\n  ", format(alt, digits = 17),
    "\n  identical:", identical(qlo, alt), "\n\n")

fd <- function(f, x, h) vapply(seq_along(x), function(j) {
  xp <- x; xp[[j]] <- xp[[j]] + h
  xm <- x; xm[[j]] <- xm[[j]] - h
  (f(xp) - f(xm)) / (2 * h)
}, 0)

# ke, k31 and ka as tape inputs; k13 is held fixed, because the
# derivative with respect to the parameter that SPLITS a double root is
# genuinely infinite there and a NaN is the honest answer for it
f <- function(k31) function(th) sum(frm_lincmt(
  parms = list(ke = exp(th[1]), k12 = K12, k21 = K21, k13 = K13,
               k31 = exp(th[2]), ka = exp(th[3]), V = 10),
  times = tt, ncmt = 3, depot = TRUE, events = ev))

cat("driving k31 onto the slow root:\n")
cat(sprintf("%12s %12s %8s %12s %12s\n", "k31/slow - 1", "eigen split",
            "finite", "rel vs FD", "FD spread"))
for (dk in c(0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)) {
  k31 <- qlo * (1 + dk)
  d <- lincmt_disp(3L, KE, K12, K21, K13, k31)
  lam <- sort(unlist(d[["lam"]]))
  th <- c(log(KE), log(k31), log(1.1))
  g <- as.numeric(RTMB::MakeTape(f(k31), th)$jacobian(th))
  g1 <- fd(f(k31), th, 1e-6)
  g2 <- fd(f(k31), th, 1e-5)
  cat(sprintf("%12.0e %12.3e %8s %12.3e %12.3e\n", dk, min(diff(lam)),
              all(is.finite(g)), max(abs(g - g1)) / max(abs(g1)),
              max(abs(g1 - g2)) / max(abs(g1))))
}

cat("\nthe VALUE either side of the line, against frm_ode() at 1e-12:\n")
d3 <- function(t, y, p) list(c(-p[6] * y[1],
  p[6] * y[1] - (p[1] + p[2] + p[4]) * y[2] + p[3] * y[3] + p[5] * y[4],
  p[2] * y[2] - p[3] * y[3], p[4] * y[2] - p[5] * y[4]))
for (k31 in c(qlo, qlo * (1 + 1e-8))) {
  a <- frm_lincmt(parms = list(ke = KE, k12 = K12, k21 = K21,
                               k13 = K13, k31 = k31, ka = 1.1,
                               V = 10),
                  times = tt, ncmt = 3, depot = TRUE, events = ev)
  o <- frm_ode(d3, init = list(0, 0, 0, 0), times = tt,
               parms = list(KE, K12, K21, K13, k31, 1.1),
               states = c("depot", "central", "peripheral1",
                          "peripheral2"),
               output = "central", events = ev, atol = 1e-12,
               rtol = 1e-12) / 10
  cat(sprintf("  k31 = %.17g\n    max difference / scale %.3e\n", k31,
              max(abs(a - o)) / max(abs(o))))
}
