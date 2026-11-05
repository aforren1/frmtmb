# Fix round 1, item 5: whole-fit time per design, one arm per process; the driver
# interleaves arms. Includes a fixed-arithmetic control so that load can be read.
# Usage: REV_ARM=after|before Rscript dev/review-impsamp-fr1-time.R
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- Sys.getenv("REV_ARM", "after")
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
ctrl <- function() {
  t0 <- proc.time()[["elapsed"]]
  s <- 0
  for (i in 1:2e7) s <- s + i * 1e-9
  proc.time()[["elapsed"]] - t0
}
for (nm in c("smallpois", "offpois", "probit", "nb", "unbal_bern", "cloglog", "wtbin", "re_sigma")) {
  de <- design_list[[nm]](1L)
  c0 <- ctrl()
  t1 <- system.time(f <- suppressWarnings(tryCatch(
    frm(de$formula, family = de$family, data = de$data, importance = 200L),
    error = function(e) e)))[["elapsed"]]
  cat(sprintf("TIME\t%s\t%s\t%.3f\t%s\t%.3f\n", arm, nm, t1,
              if (inherits(f, "error")) "stop" else f$importance$rounds, c0))
}
