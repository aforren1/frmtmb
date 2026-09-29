## Reviewer, claim 8: what a refit loop pays for the new stage.
##
## The worker timed the two BUILDS in two processes. This times the two
## readings inside ONE process on the lane build, so the arms share
## every other code path: arm "tol" is the default grad_tol, arm "off"
## sets grad_tol above the fit's own gradient so no Hessian is built.
## restarts = 0 in both, so grad_tol governs only the warning path. Arms
## are interleaved per round, the figure is the minimum over rounds, and
## the control is a second "off" arm that must read 1.00.
## usage: Rscript gradcheck-rev-14-boot.R <core-lib> [rounds] [nsim]
a <- commandArgs(TRUE)
LIB <- a[1]
ROUNDS <- if (length(a) >= 2) as.integer(a[2]) else 5L
NSIM <- if (length(a) >= 3) as.integer(a[3]) else 20L
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), " rounds:", ROUNDS,
    " nsim:", NSIM, "\n\n")
`%||%` <- function(x, y) if (is.null(x)) y else x

run <- function(lbl, mk) {
  f0 <- suppressWarnings(mk(1e-3))
  g0 <- max(abs(drop(f0$obj$gr(f0$opt$par)) * (f0$par_units %||% 1)))
  off <- max(10 * g0, 1)
  tt <- to <- tc <- rep(NA_real_, ROUNDS)
  bt <- function(gt) {
    f <- suppressWarnings(mk(gt))
    t0 <- proc.time()[["elapsed"]]
    invisible(suppressWarnings(frm_bootstrap(f, nsim = NSIM, seed = 1L)))
    proc.time()[["elapsed"]] - t0
  }
  for (k in seq_len(ROUNDS)) {
    tt[k] <- bt(1e-3); to[k] <- bt(off); tc[k] <- bt(off)
  }
  cat(sprintf(paste0("%-32s np %2d gmax %9.3e trips %-5s  ",
                     "tol %7.3f s  off %7.3f s  ratio %5.3f  ",
                     "CONTROL %5.3f\n"),
              lbl, length(f0$opt$par), g0, g0 > 1e-3,
              min(tt), min(to), min(tt) / min(to), min(tc) / min(to)))
  cat(sprintf("%-32s per-round tol %s\n", "",
              paste(format(tt, digits = 4), collapse = " ")))
  cat(sprintf("%-32s per-round off %s\n", "",
              paste(format(to, digits = 4), collapse = " ")))
  cat(sprintf("%-32s per-round ctl %s\n\n", "",
              paste(format(tc, digits = 4), collapse = " ")))
}

ctl <- function(gt) frmtmb_control(restarts = 0, grad_tol = gt)

set.seed(1040)
n <- 20000
d1 <- data.frame(x = rnorm(n), z = rnorm(n))
d1$y <- rnorm(n, 1 + 2 * d1$x - 0.5 * d1$z, 2)
run("gaussian n=20000 (trips)",
    function(gt) frm(bf(y ~ x + z), family = gaussian(), data = d1,
                     control = ctl(gt)))

set.seed(1041)
n <- 400
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$y <- rnorm(n, 1 + 2 * d2$x - 0.5 * d2$z, 2)
run("gaussian n=400 (does not trip)",
    function(gt) frm(bf(y ~ x + z), family = gaussian(), data = d2,
                     control = ctl(gt)))

set.seed(1042)
n <- 20000
d3 <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
e3 <- 0.8 * d3$x1 - 0.5 * d3$x2
d3$yo <- cut(e3 + rlogis(n), breaks = c(-Inf, -1, 0.5, 2, Inf),
             labels = FALSE)
run("cumulative n=20000 (trips)",
    function(gt) frm(bf(yo ~ x1 + x2), family = cumulative(), data = d3,
                     control = ctl(gt)))
cat("DONE boot\n")
