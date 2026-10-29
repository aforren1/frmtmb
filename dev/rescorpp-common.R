## Lane rescorpp helpers. First argument: "lane" (default) or "base".
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
base <- c("C:/Users/adf44/source/r/rellib-r15",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-rescorpp-lib", base)
} else base)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/stan-cache")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("arm:", arm, " lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "| sample",
    dirname(find.package("frmtmb.sample")),
    format(packageVersion("frmtmb.sample")), "\n")
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
