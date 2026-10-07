# Lane optima, nlminb's reported point: on how many fits of
# dev/optima-csmix.R's mixture designs is the objective AT the reported
# estimates not the objective the fit reports (logLik)? On 0.68.1 nlminb
# could return the rejected trial as `par` beside the best value.
#   Rscript dev/optima-parcheck.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
designs <- list(
  mix_cum_sratio = function(d) frm(bf(y ~ x + cs(z)),
                                   family = mixture(cumulative(), sratio()),
                                   data = d),
  mix_cum_cum = function(d) frm(bf(y ~ x + cs(z)),
                                family = mixture(cumulative(), cumulative()),
                                data = d),
  mix_sratio_sratio = function(d) frm(bf(y ~ x + cs(z)),
                                      family = mixture(sratio(), sratio()),
                                      data = d),
  mix_cum_sratio_nocs = function(d) frm(bf(y ~ x + z),
                                        family = mixture(cumulative(),
                                                         sratio()),
                                        data = d))
for (nm in names(designs)) {
  res <- vapply(1:20, function(s) {
    f <- tryCatch(suppressWarnings(designs[[nm]](mk(s, mix = TRUE))),
                  error = function(e) NULL)
    if (is.null(f)) return("error")
    v <- f$obj$fn(f$opt$par)
    if (!is.finite(v)) "par not finite" else
      if (v != f$opt$objective) {
        cat(sprintf(paste0("  %s seed %d: objective at par - reported = ",
                           "%.3g (code %d, %s)\n"),
                    nm, s, v - f$opt$objective, f$opt$convergence,
                    f$opt$message))
        "par elsewhere"
      } else "par at objective"
  }, "")
  tb <- table(res)
  cat(sprintf("%-22s %s\n", nm, paste(names(tb), tb, sep = "=",
                                      collapse = " ")))
}
