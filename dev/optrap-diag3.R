# Lane optrap, item 1: what code 1 at y x 1e6 is (dev/optrap-scale.R,
# design "slope", seeds 2 and 22; "ri20" seed 37).
#   Rscript dev/optrap-diag3.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(if (nzchar(args[1])) args[1],
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", find.package("frmtmb"), "\n")
src <- readLines("dev/optrap-scale.R")
eval(parse(text = src[grep("^mk <- function", src):
                        (grep("^X <- list", src) - 1L)]))
for (cs in list(list("slope", 2), list("slope", 22), list("ri20", 37))) {
  m <- mk(cs[[1]], cs[[2]])
  for (sc in c(1, 1e6)) {
    d <- m$d
    d$y <- d$y * sc
    w <- character()
    f <- withCallingHandlers(frm(m$f, data = d, verbose = 2),
      warning = function(x) {w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")},
      message = function(x) {
        if (!grepl("^ *[0-9]+:", conditionMessage(x))) {
          cat("  M:", substr(conditionMessage(x), 1, 150), "\n")
        }
        invokeRestart("muffleMessage")})
    cat(sprintf("%s seed %d x %g: code %d (%s) obj %.10g iter %s eval %s\n",
                cs[[1]], cs[[2]], sc, f$opt$convergence, f$opt$message,
                f$opt$objective, f$opt$iterations,
                paste(f$opt$evaluations, collapse = "/")))
    cat("  par/units:", signif(f$opt$par / (f$par_units %||% 1), 6), "\n")
    cat("  grad*units:", signif(f$obj$gr(f$opt$par) * (f$par_units %||% 1),
                                3), "\n")
    for (x in w) cat("   W:", substr(x, 1, 120), "\n")
  }
}
