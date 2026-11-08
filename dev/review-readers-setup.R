# Library setup for the adversarial review of lane readers.
# REV_ARM=before reads rellib-r18 alone; REV_ARM=mutant puts rev-readers-lib
# first; anything else puts the lane library first.
arm <- Sys.getenv("REV_ARM", "after")
lib_lane <- "C:/Users/adf44/source/r/wt-readers-lib"
lib_base <- "C:/Users/adf44/source/r/rellib-r18"
lib_user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
lib_mut <- "C:/Users/adf44/source/r/rev-readers-lib"
.libPaths(switch(arm,
  before = c(lib_base, lib_user),
  mutant = c(lib_mut, lib_base, lib_user),
  c(lib_lane, lib_base, lib_user)))
suppressPackageStartupMessages(library(frmtmb))
cat("arm:", arm, " lib:", dirname(find.package("frmtmb")),
    " frmtmb", format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
# count warnings by class while returning the value
with_warn <- function(expr) {
  ws <- list()
  val <- withCallingHandlers(expr, warning = function(w) {
    ws[[length(ws) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws)
}
wclass <- function(ws) {
  vapply(ws, function(w) paste(setdiff(class(w), c("warning", "condition",
                                                   "simpleWarning")),
                               collapse = "/"), "")
}
