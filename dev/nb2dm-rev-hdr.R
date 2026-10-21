# Reviewer header: the reviewer's library (lane build), rellib-r10, user.
# Source it first. NB2DM_REV_ARM=base drops the reviewer's library.
arm <- Sys.getenv("NB2DM_REV_ARM", "lane")
.libPaths(unique(c(
  if (!identical(arm, "base")) "C:/Users/adf44/source/r/wt-nb2dm-rev-lib",
  "C:/Users/adf44/source/r/rellib-r10",
  "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/stan-cache",
           NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")),
    format(utils::packageVersion("frmtmb")), "arm", arm, "\n")
# run one probe, print its value or its condition, never stop
probe <- function(label, expr) {
  cat("\n### ", label, "\n", sep = "")
  ws <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      cat("ERROR: ", conditionMessage(e), "\n", sep = "")
      invisible(NULL)
    }),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cat("MESSAGE: ", conditionMessage(m), sep = "")
      invokeRestart("muffleMessage")
    })
  if (length(ws)) cat("WARNING: ", ws, sep = "\n  ")
  if (!is.null(r)) print(r)
  invisible(r)
}
