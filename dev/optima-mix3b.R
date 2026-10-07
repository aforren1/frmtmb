# Lane optima, item 2: the review's three-component ordinal mixture with
# sratio's link changed from cloglog to logit, on the arm given: does
# the fit end with a finite gradient, and what does it warn?
#   Rscript dev/optima-mix3b.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
src <- readLines("dev/optima-mix3.R")
eval(parse(text = src[grep("^set.seed", src):grep("^d <- data.frame", src)]))
for (fam in list(mixture(cumulative("probit"), sratio(), acat()),
                 mixture(cumulative("probit"), acat()),
                 mixture(cumulative("probit"), cumulative("probit")))) {
  w <- character()
  r <- tryCatch(withCallingHandlers(
    frm(bf(y ~ x), family = fam, data = d),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 90))
      invokeRestart("muffleWarning")
    }), error = function(e) conditionMessage(e))
  cat("==", fam$family, paste(fam$mix_names %||% "", collapse = ","), "\n")
  if (is.character(r)) cat("ERROR:", substr(r, 1, 100), "\n") else {
    cat("logLik", format(as.numeric(logLik(r)), digits = 10), "code",
        r$opt$convergence, "gradient finite",
        all(is.finite(r$obj$gr(r$opt$par))), "\n")
  }
  if (length(w)) cat(paste("  warning:", w), sep = "\n")
}
