## Reviewer, claim 1(b) and claim 8: what stage 2 costs on a fit that
## trips, measured WITHIN one build so the two arms share every other
## code path.
##
## Instrument. Both arms set restarts = 0, so grad_tol governs only the
## warning path and not the optimizer. Arm "tol" is the default
## grad_tol; arm "off" sets grad_tol above the fit's own gradient, so
## the trip-wire never trips and no Hessian is built. The arms are
## interleaved inside one round and the reported figure is the MINIMUM
## over rounds. The control is a second "off" arm, which must read 1.00.
## usage: Rscript gradcheck-rev-05-cost.R <core-lib> [rounds]
a <- commandArgs(TRUE)
LIB <- a[1]
ROUNDS <- if (length(a) >= 2) as.integer(a[2]) else 3L
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "  rounds:", ROUNDS, "\n")
cat("proc.time() resolution probe: ")
t0 <- proc.time()[["elapsed"]]
repeat { t1 <- proc.time()[["elapsed"]]; if (t1 > t0) break }
cat(format((t1 - t0) * 1000, digits = 3), "ms\n\n")

`%||%` <- function(x, y) if (is.null(x)) y else x
quiet <- function(expr) suppressWarnings(expr)
ctl <- function(gt) frmtmb_control(restarts = 0, grad_tol = gt)

make <- list()
make$glmm5000 <- function() {
  set.seed(8001)
  ng <- 5000
  dd <- data.frame(g = factor(rep(seq_len(ng), 4)))
  dd$x <- rnorm(nrow(dd))
  re <- rnorm(ng, 0, 0.7)
  dd$y <- rpois(nrow(dd), exp(0.3 + 0.5 * dd$x + re[dd$g]))
  dd
}
make$glmm5000_slope <- function() {
  set.seed(8002)
  ng <- 5000
  dd <- data.frame(g = factor(rep(seq_len(ng), 4)))
  dd$x <- rnorm(nrow(dd))
  re <- rnorm(ng, 0, 0.7); rs <- rnorm(ng, 0, 0.4)
  dd$y <- rpois(nrow(dd), exp(0.3 + (0.5 + rs[dd$g]) * dd$x + re[dd$g]))
  dd
}
make$gp <- function() {
  set.seed(8003)
  n <- 600
  dd <- data.frame(t = sort(runif(n, 0, 12)))
  dd$y <- rnorm(n, sin(dd$t) + 0.2 * dd$t, 0.3)
  dd
}
make$smooth <- function() {
  set.seed(8004)
  n <- 20000
  dd <- data.frame(x = sort(runif(n, 0, 8)), z = rnorm(20000))
  dd$y <- rnorm(n, sin(dd$x) + 0.4 * dd$z, 0.5)
  dd
}
fitf <- list(
  glmm5000 = function(dd, gt) quiet(frm(bf(y ~ x + (1 | g)),
                                        family = poisson(), data = dd,
                                        control = ctl(gt))),
  glmm5000_slope = function(dd, gt) quiet(frm(bf(y ~ x + (x | g)),
                                              family = poisson(), data = dd,
                                              control = ctl(gt))),
  gp = function(dd, gt) quiet(frm(bf(y ~ gp(t, k = 30)),
                                  family = gaussian(),
                                  data = dd, control = ctl(gt))),
  smooth = function(dd, gt) quiet(frm(bf(y ~ s(x) + z), family = gaussian(),
                                      data = dd, control = ctl(gt))))

tm <- function(f) { t0 <- proc.time()[["elapsed"]]
  v <- f(); list(t = proc.time()[["elapsed"]] - t0, v = v) }

for (nm in names(make)) {
  dd <- make[[nm]]()
  ff <- fitf[[nm]]
  # one throwaway fit to learn the gradient, so the "off" arm's grad_tol
  # is above it by construction rather than by guess
  f0 <- ff(dd, 1e-3)
  g0 <- max(abs(drop(f0$obj$gr(f0$opt$par)) * (f0$par_units %||% 1)))
  np <- length(f0$opt$par)
  d0 <- diagnose(f0, quiet = TRUE)
  off_tol <- max(10 * g0, 1)
  trips <- g0 > 1e-3
  ttol <- toff <- tctl <- rep(NA_real_, ROUNDS)
  for (k in seq_len(ROUNDS)) {
    ttol[k] <- tm(function() ff(dd, 1e-3))$t
    toff[k] <- tm(function() ff(dd, off_tol))$t
    tctl[k] <- tm(function() ff(dd, off_tol))$t
  }
  cat(sprintf(
    "%-16s np=%-3d n=%-6d gmax=%-10s trips=%-5s headroom=%-10s\n",
    nm, np, nrow(dd), format(g0, digits = 4), trips,
    format(d0$grad_headroom %||% NA, digits = 4)))
  cat(sprintf(
    "%-16s fit+stage2 %.3f s   fit only %.3f s   ratio %.3f   CONTROL %.3f\n",
    "", min(ttol), min(toff), min(ttol) / min(toff),
    min(tctl) / min(toff)))
  cat(sprintf("%-16s per-round tol: %s\n", "",
              paste(format(ttol, digits = 4), collapse = " ")))
  cat(sprintf("%-16s per-round off: %s\n", "",
              paste(format(toff, digits = 4), collapse = " ")))
  cat(sprintf("%-16s per-round ctl: %s\n\n", "",
              paste(format(tctl, digits = 4), collapse = " ")))
}
cat("DONE cost\n")
