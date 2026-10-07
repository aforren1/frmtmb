# Lane optima, item 1: wall-clock cost of a mo() fit, base against lane,
# with a control. One arm per process; each process times the 50 fits
# of seeds 1..50 (ls ~ mo(income) * age, and ls ~ mo(income)) and a
# fixed control computation built from the same R code on both arms
# (QR solves of a fixed matrix), so a loaded machine shows up in the
# control rather than in the ratio.
#   Rscript dev/optima-mo-time.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
src <- readLines("dev/optima-mo-study.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^# exact maximum", src) - 1L)]))
ds <- lapply(1:50, mk)
invisible(suppressWarnings(frm(ls ~ mo(income) * age, data = ds[[1]])))
control <- function() {
  set.seed(1)
  M <- matrix(rnorm(250000), 500)
  t0 <- proc.time()[[3]]
  for (i in 1:40) qr.solve(M, M[, 1])
  proc.time()[[3]] - t0
}
tm <- function(fo) {
  t0 <- proc.time()[[3]]
  for (d in ds) suppressWarnings(frm(fo, data = d))
  proc.time()[[3]] - t0
}
c1 <- control()
ti <- tm(ls ~ mo(income) * age)
tmn <- tm(ls ~ mo(income))
c2 <- control()
cat(sprintf("TIME arm %s interaction %.2f s main %.2f s control %.2f %.2f s\n",
            arm, ti, tmn, c1, c2))
