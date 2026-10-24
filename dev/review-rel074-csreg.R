# Release review of 0.74.0: frmtmb_register_covstruct()'s name
# refusals, one name per call, printed before each call.
#   Rscript dev/review-rel074-csreg.R [names...]
.libPaths(c("C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
d <- function(...) stop("x")
a <- commandArgs(trailingOnly = TRUE)
nms <- if (length(a)) a else c("gaussian", "dec", "poly", "ns", "log",
                               "car", "mystruct")
for (nm in nms) {
  cat("calling", nm, "\n"); flush.console()
  r <- tryCatch({
    frmtmb_register_covstruct(nm, lpdf = d, vcov = d,
                              npar = function(...) 1L,
                              package = "frmtmb.test")
    "registered"
  }, error = function(e) substr(conditionMessage(e), 1, 220))
  cat(nm, ":", r, "\n"); flush.console()
}
