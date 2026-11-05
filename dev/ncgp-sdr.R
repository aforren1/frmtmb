# Which part of sdreport() costs what on an exact gp()? Seed 20261103.
#   Rscript dev/ncgp-sdr.R <n> <what>
a <- commandArgs(TRUE)
n <- as.integer(a[1]); what <- a[2]
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261103)
x <- sort(runif(n, 0, 10))
d <- data.frame(x = x, y = sin(x) + 0.3 * rnorm(n))
fr <- frmtmb:::frm(bf(y ~ gp(x)), data = d, family = gaussian(), dry_run = "frame")
tpl <- fr$par_template
pk <- function(tag) {
  p <- ps <- system2("powershell", c("-NoProfile", "-c", sprintf("(Get-Process -Id %d).PeakWorkingSet64/1MB", Sys.getpid())), stdout = TRUE)
  cat(sprintf("NCGPSR n=%d %s peakMB %s t=%.2f\n", n, tag, p, proc.time()[[3]] - t0))
}
t0 <- proc.time()[[3]]
obj <- frmtmb:::gp_tape_scope(fr, RTMB::MakeADFun(frmtmb:::build_objective(fr), tpl, random = "b", silent = TRUE))
pk("tape")
op <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 100))
pk("opt")
scope <- function(e) if (what == "scope") frmtmb:::gp_tape_scope(fr, e) else if (what == "atomic") { o <- RTMB::TapeConfig(); RTMB::TapeConfig(matmul = "atomic"); on.exit(RTMB::TapeConfig(o)); e } else e
s1 <- scope(RTMB::sdreport(obj, par.fixed = op$par, getJointPrecision = FALSE, skip.delta.method = TRUE))
pk("sdreport no joint")
s2 <- scope(RTMB::sdreport(obj, par.fixed = op$par, getJointPrecision = TRUE, skip.delta.method = TRUE))
pk("sdreport joint")
print(names(s2))
print(dim(s2$jointPrecision))
