# Reviewer 2: does evaluating lcdf/lccdf on every row (lib7) cost more
# tape memory than on the censored rows only (lane)? Peak working set of
# this process after taping and one fn + gr evaluation, on a 10 x 300
# cleft design (the worker's cleft arm, scaled down). No optimization.
# Usage: Rscript r2-peakmem.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[[1]], "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
pk <- function() as.numeric(system(sprintf(
  "powershell -NoProfile -Command \"(Get-Process -Id %d).PeakWorkingSet64\"", Sys.getpid()),
  intern = TRUE)) / 2^20
m0 <- pk()
set.seed(1); NS <- 10L; NT <- 300L
u <- rnorm(NS, 0, 0.35); b <- rnorm(NS, 0, 0.2)
d <- do.call(rbind, lapply(seq_len(NS), function(s) {
  cond <- rep(0:1, length.out = NT)
  x <- ddm_simulate(NT, mu = 0.4 + 0.9 * cond + u[s], bs = 1.4 * exp(b[s]), ndt = 0.25, bias = 0.5)
  x$cond <- cond; x$s <- factor(s); x }))
d$code <- 0L; d$y2 <- d$rt
fast <- d$rt < 0.45; d$code[fast] <- -1L; d$rt[fast] <- 0.45
pick <- which(!fast)[seq(5, sum(!fast), by = 5)]
lo <- floor(d$rt[pick] * 10) / 10
d$code[pick] <- 2L; d$y2[pick] <- lo + 0.1; d$rt[pick] <- pmax(lo, 0.45)
t0 <- proc.time()[["elapsed"]]
o <- frm(bf(rt | dec(upper) + cens(code, y2) ~ cond + (1 | s), bs ~ 1 + (1 | s), ndt ~ 1,
            bias = 0.5), family = wiener(), data = d, dry_run = "objective")
f <- o$obj$fn(o$obj$par); g <- o$obj$gr(o$obj$par)
cat(sprintf("%s: rows %d (left %d, interval %d); fn %.10f; peak MB before %.0f after %.0f; %.1f s\n",
            basename(dirname(find.package("frmtmb.eam"))), nrow(d), sum(d$code == -1),
            sum(d$code == 2), f, m0, pk(), proc.time()[["elapsed"]] - t0))
