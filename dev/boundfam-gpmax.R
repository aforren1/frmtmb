# Lane boundfam, consolidation item 2: is 0.81.0's point of the two
# student rows a proper maximum? Refit from the base estimates on the lane
# build and print the conditions, the stop and the standard errors.
#   Rscript dev/boundfam-gpmax.R <lane lib> <base .rds of dev/boundfam-gp.R>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-rel082/dev/gpstop-fx2-sweep.R")
eval(parse(text = src[45:87]))
b <- readRDS(a[2])
for (nm in names(b)) {
  body <- sub("_n.*", "", nm)
  d <- b[[nm]]$data
  f0 <- suppressWarnings(suppressMessages(frm(spell(body), data = d,
                                               family = student())))
  tpl <- f0$frame[["par_template"]]
  p <- b[[nm]]$par
  st <- list()
  for (cp in unique(names(p))) st[[cp]] <- unname(p[names(p) == cp])
  w <- character(); m <- character()
  f <- withCallingHandlers(frm(spell(body), data = d, family = student(),
                               start = st),
    warning = function(e) { w <<- c(w, substr(conditionMessage(e), 1, 120)); invokeRestart("muffleWarning") },
    message = function(e) { m <<- c(m, substr(conditionMessage(e), 1, 120)); invokeRestart("muffleMessage") })
  cat(sprintf("\n== %s from 0.81.0's point: ll %.6f (base %.6f) code %d max|grad| %.3g\n",
              nm, as.numeric(logLik(f)), b[[nm]]$ll, f$opt$convergence,
              max(abs(f$obj$gr(f$opt$par)))))
  for (x in w) cat("  W:", x, "\n")
  for (x in m) cat("  M:", x, "\n")
  se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
  cat("  SE:", sprintf("%.4g", se), "\n")
}
