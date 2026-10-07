# Lane optima, item 3: where in frm()'s own optimizer path does the
# "NA/NaN gradient evaluation" of a cs() ordinal mixture arise? Wraps
# the taped objective's gr() to record every gradient call (point,
# objective there, finite or not) during a real frm() fit.
#   Rscript dev/optima-csmix-where.R base|lane seed design
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
seed <- as.integer(args[2])
design <- args[3]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
d <- mk(seed, mix = TRUE)
fam <- switch(design,
  cum_sratio = mixture(cumulative(), sratio()),
  cum_cum = mixture(cumulative(), cumulative()))
log <- list()
trace(RTMB::MakeADFun, exit = quote({
  o <- returnValue()
  g0 <- o$gr
  f0 <- o$fn
  o$gr <- function(x, ...) {
    g <- g0(x, ...)
    f <- f0(x)
    log[[length(log) + 1L]] <<- list(x = x, f = f, ok = all(is.finite(g)),
                                     g = g)
    g
  }
  assign("o_last", o, envir = globalenv())
}), print = FALSE)
r <- tryCatch(suppressWarnings(frm(bf(y ~ x + cs(z)), family = fam, data = d,
                                   verbose = 1)),
              error = function(e) conditionMessage(e))
untrace(RTMB::MakeADFun)
if (is.character(r)) cat("ERROR:", r, "\n")
ok <- vapply(log, `[[`, NA, "ok")
fv <- vapply(log, `[[`, 0, "f")
cat("gradient calls", length(log), "; non-finite gradients", sum(!ok),
    "; at a finite objective", sum(!ok & is.finite(fv)), "\n")
if (any(!ok)) {
  b <- log[[which(!ok)[1]]]
  cat("first: objective", format(b$f, digits = 12), "\n")
  print(b$g)
  print(b$x)
  saveRDS(b, sprintf("dev/optima-log/csmix-where-%s-%d.rds", design, seed))
}
