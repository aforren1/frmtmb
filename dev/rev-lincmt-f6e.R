# Round 2, item 1, part five, and the question everything turns on.
#
# dev/rev-lincmt-f6d.R found two different regimes behind the same
# NaN:
#
#   1  a GENUINE double root at ordinary scale (k13 = 1e-300, k31 on
#      the slow root). The value is right, the true gradient is finite,
#      and clamping `arg` recovers it at no cost to the value.
#   2  a wide SPREAD of rate constants, eight decades or more, where
#      `arg` rounds to one because `pp` and `qq` have lost their
#      relative precision. There the clamp changes the value by a
#      factor of 1e2 to 1e6, which says the two small eigenvalues were
#      already unresolvable.
#
# In regime 2 the only question that matters is whether the VALUE the
# shipped code returns is right. If it is, the NaN gradient is a loud
# refusal protecting a right answer and leaving it is defensible. If it
# is not, the NaN is the least of the problem.
#
# Script path: dev/rev-lincmt-f6e.R. Seed 4242, the same draws
# dev/rev-lincmt-f6d.R used.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-mpfr.R")

tt <- c(0.5, 1, 2, 4, 8, 12, 24, 48)
ev <- data.frame(time = 0, state = "depot", value = 100, ii = 8,
                 addl = 3L)
gradof <- function(p) {
  f <- function(th) sum(frm_lincmt(
    parms = list(ke = exp(th[1]), k12 = exp(th[2]), k21 = exp(th[3]),
                 k13 = exp(th[4]), k31 = exp(th[5]), ka = exp(th[6]),
                 V = 10),
    times = tt, ncmt = 3, depot = TRUE, events = ev))
  x <- log(p)
  as.numeric(MakeTape(f, x)$jacobian(x))
}
# one bolus into the depot, no schedule, so the reference is one
# matrix exponential per lag
lin1 <- function(p, u) frm_lincmt(
  parms = list(ke = p[1], k12 = p[2], k21 = p[3], k13 = p[4],
               k31 = p[5], ka = p[6], V = 1),
  times = u, ncmt = 3, depot = TRUE, init = list(depot = 1),
  output = "central")

set.seed(4242)
bad <- list()
while (length(bad) < 25L) {
  p <- 0.2 * 10^runif(6, -5, 5)
  if (!all(is.finite(gradof(p)))) bad[[length(bad) + 1L]] <- p
}
cat("\n=== the VALUE where the gradient is NaN ===\n")
cat("25 draws from the 10-decade box that make the tape's gradient",
    "\nNaN. The value at seven lags against the 300-bit reference,",
    "\nrelative to the trajectory's own maximum.\n\n")
lags <- c(0.05, 0.5, 2, 8, 24, 72, 200)
worst <- 0; wat <- NULL; n <- 0L
errs <- numeric(0)
for (p in bad) {
  a <- lin1(p, lags)
  b <- vapply(lags, function(u) as.numeric(
    bolus_ss(3L, TRUE, list(ke = p[1], k12 = p[2], k21 = p[3],
                            k13 = p[4], k31 = p[5], ka = p[6]), u)), 0)
  pk <- max(abs(b))
  if (!(pk > 0) || anyNA(b)) next
  e <- max(abs(a - b)) / pk
  n <- n + 1L
  errs <- c(errs, e)
  if (e > worst) { worst <- e; wat <- p }
}
cat("  cases compared:", n, "\n")
cat("  median error / peak:", format(median(errs), digits = 4), "\n")
cat("  worst  error / peak:", format(worst, digits = 4), "\n")
cat("  at ", paste(format(wat, digits = 4), collapse = "  "), "\n")
cat("\n  the distribution:\n")
print(signif(quantile(errs, c(0, 0.25, 0.5, 0.75, 0.9, 1)), 4))

cat("\n=== the same draws, but the gradient FINITE, as a control ===\n")
set.seed(4242)
good <- list()
while (length(good) < 25L) {
  p <- 0.2 * 10^runif(6, -5, 5)
  if (all(is.finite(gradof(p)))) good[[length(good) + 1L]] <- p
}
g <- numeric(0)
for (p in good) {
  a <- lin1(p, lags)
  b <- vapply(lags, function(u) as.numeric(
    bolus_ss(3L, TRUE, list(ke = p[1], k12 = p[2], k21 = p[3],
                            k13 = p[4], k31 = p[5], ka = p[6]), u)), 0)
  pk <- max(abs(b))
  if (!(pk > 0) || anyNA(b)) next
  g <- c(g, max(abs(a - b)) / pk)
}
cat("  cases compared:", length(g), "  worst error / peak:",
    format(max(g), digits = 4), "\n")
cat("\nIf the two distributions are the same, the NaN gradient is",
    "\nguarding a value that is right and the refusal is protective.",
    "\nIf the NaN set is worse, the NaN is a symptom.\n")
