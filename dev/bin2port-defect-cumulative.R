# The cumulative tied-threshold defect found by the value check of
# log_lik:462 (lane bin2port). frmtmb 0.71.0's cumulative log
# probability of an END category is NaN when the two thresholds of the
# interior pair it is clamped to are equal, though that category has
# positive probability. Log: dev/bin2port-out/defect-cumulative.txt.
#
#   Rscript dev/bin2port-defect-cumulative.R [lib]
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/rellib-r9"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
lp <- frmtmb::cumulative()[["ord_lpdf_make"]](identity)
dens <- function(y, mu, tau) {
  as.numeric(lp(y, list(mu = mu, disc = 1), list(), list(tau_raw = tau)))
}
# the exact log probabilities under the logit link
exact <- function(y, mu, tau) {
  vapply(seq_along(mu), function(s) {
    F <- c(0, stats::plogis(tau - mu[s]), 1)
    log(F[y[s] + 1L] - F[y[s]])
  }, 0)
}
mu <- c(-0.7, 0, 1.3)
for (tau in list(c(0, 0, 0), c(-1, -1, 1), c(-1, 1, 1), c(-1, 0, 1))) {
  cat("thresholds", format(tau), "\n")
  for (y in 1:4) {
    a <- dens(rep(y, 3), mu, tau)
    e <- exact(rep(y, 3), mu, tau)
    cat(sprintf("  y = %d  frmtmb %s   exact %s\n", y,
                paste(format(a, digits = 6), collapse = " "),
                paste(format(e, digits = 6), collapse = " ")))
  }
}
# The same on the RTMB tape, through a fit's objective: a cumulative
# model whose two lowest thresholds are held equal by a map is not
# something frm() builds, so the tape is checked on the density alone:
# taped at tied thresholds, is category 1 NaN there too?
f <- function(p) {
  tau <- c(p[1], p[1], p[2])
  sum(lp(c(1, 1, 4), list(mu = c(-0.7, 0, 1.3), disc = 1), list(),
         list(tau_raw = tau)))
}
tp <- RTMB::MakeTape(f, c(0, 1))
cat("taped sum over rows y = 1, 1, 4 at tau = (0, 0, 1):", tp(c(0, 1)),
    "\n")
cat("exact:", sum(exact(c(1, 1, 4), c(-0.7, 0, 1.3), c(0, 0, 1))), "\n")
cat("on doubles:", f(c(0, 1)), "\n")

# hurdle_cumulative (review 2026-10-14, finding 4): it calls the same
# function on pmax(y, 1) and masks it with (1 - i0), so a y = 0 row is
# NaN too at a tie of the lowest threshold pair
hp <- frmtmb::hurdle_cumulative()[["ord_lpdf_make"]](identity)
hdens <- function(y, mu, hu, tau) {
  as.numeric(hp(y, list(mu = mu, hu = hu, disc = 1), list(),
                list(tau_raw = tau)))
}
hexact <- function(y, mu, hu, tau) {
  vapply(seq_along(y), function(s) {
    if (y[s] == 0) return(log(hu))
    log(1 - hu) + exact(y[s], mu[s], tau)
  }, 0)
}
y <- c(0, 1, 4)
mu <- c(0, 0, 0)
for (tau in list(c(0, 0, 1), c(-1, 0, 1))) {
  cat("hurdle_cumulative, hu 0.2, y", y, ", thresholds", format(tau), "\n")
  cat("  frmtmb", format(hdens(y, mu, 0.2, tau), digits = 7), "\n")
  cat("  exact ", format(hexact(y, mu, 0.2, tau), digits = 7), "\n")
}
