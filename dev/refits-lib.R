# Shared setup for the refits lane's scripts. ARM=before runs on the base
# build (rellib-r21) alone; anything else puts the lane library first.
arm <- Sys.getenv("ARM", "after")
lanelib <- "C:/Users/adf44/source/r/wt-refits-lib"
base <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "before")) base else c(lanelib, base))
suppressPackageStartupMessages(library(frmtmb))
cat("arm:", arm, " lib:", find.package("frmtmb"),
    as.character(packageVersion("frmtmb")), "\n")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)), 1L,
                            500L)
go <- function(lab, expr) {
  cat("\n-- ", lab, " --\n", sep = "")
  ws <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e)
      paste0("ERROR [", class(e)[1L], "]: ", short(e))),
    warning = function(w) {
      ws <<- c(ws, paste("WARN:", short(w))); invokeRestart("muffleWarning")
    },
    message = function(m) {
      ws <<- c(ws, paste("MSG:", short(m))); invokeRestart("muffleMessage")
    })
  if (length(ws)) cat("conditions:\n", paste(" *", ws, collapse = "\n"), "\n")
  print(r)
  invisible(r)
}
