# Where in core a saved fit's absent-package family is first called:
# the call stack at the error, for each failing method.
#   Rscript dev/splitsurv-saved-trace.R <dir> <scratch lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[2], "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
s <- readRDS(file.path(a[1], "saved-fits.rds"))
drop <- c("tryCatch", "withCallingHandlers", "doTryCatch", "tryCatchList",
          "tryCatchOne", "force", "trace", ".handleSimpleError", "h")
trace <- function(lab, expr) {
  calls <- NULL
  tryCatch(withCallingHandlers(force(expr), error = function(e) {
    calls <<- vapply(sys.calls(), function(cl) {
      paste(deparse(cl[[1]])[1], collapse = "")
    }, "")
  }), error = function(e) NULL)
  cat("==", lab, ":", paste(calls[!calls %in% drop], collapse = " > "),
      "\n")
}
trace("summary cox", summary(s$fc))
trace("vcov cox", vcov(s$fc))
trace("update cox", update(s$fc, . ~ 1))
trace("simulate rp", simulate(s$fr, nsim = 1))
trace("update rp", update(s$fr, . ~ 1))
