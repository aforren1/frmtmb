## Reviewer helpers, sourced by dev/review-mointerop-*.R.
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
base <- c("C:/Users/adf44/source/r/rellib-r14",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mointerop-lib", base)
} else base)
suppressMessages({
  library(frmtmb)
  library(marginaleffects)
  library(emmeans)
})
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "\n")
options(width = 120)
show <- function(tag, expr) {
  cat("\n==", tag, "==\n")
  w <- character(0)
  r <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                           warning = function(x) {
                             w <<- c(w, conditionMessage(x))
                             invokeRestart("muffleWarning")
                           })
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n")
  else print(r)
  cat("warnings:", length(w), "\n")
  for (m in w) cat("WARNING:", substr(m, 1, 300), "\n")
  invisible(r)
}
