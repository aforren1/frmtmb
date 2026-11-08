# Review of lane arstep: shared helpers. "lane" is wt-arstep-lib over
# rellib-r18; "base" is rellib-r18 alone (frmtmb 0.79.1 at 3cd2a0c8); "rev" is
# a mutant build in rev-arstep-lib over rellib-r18.
RV_LANE <- "C:/Users/adf44/source/r/wt-arstep-lib"
RV_BASE <- "C:/Users/adf44/source/r/rellib-r18"
RV_REV <- "C:/Users/adf44/source/r/rev-arstep-lib"
RV_USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
RV_WT <- "C:/Users/adf44/source/r/frmtmb-wt-arstep"
RV_OUT <- file.path(RV_WT, "dev/review-arstep-out")
dir.create(RV_OUT, showWarnings = FALSE)
rv_libs <- function(arm = c("lane", "base", "rev")) {
  arm <- match.arg(arm)
  .libPaths(c(switch(arm, lane = RV_LANE, rev = RV_REV, base = NULL),
              RV_BASE, RV_USER))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
             FRMTMB_STAN_CACHE = file.path(RV_WT, "dev/stan-cache"),
             NOT_CRAN = "true")
  suppressPackageStartupMessages(library(frmtmb))
  cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
      as.character(packageVersion("frmtmb")), "\n")
  invisible(arm)
}
rv_put <- function(name, f) {
  ns <- asNamespace("frmtmb")
  environment(f) <- ns
  unlockBinding(name, ns); assign(name, f, envir = ns); lockBinding(name, ns)
}
# warnings of a call, as text, and the value
rv_quiet <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  list(value = v, warnings = w)
}
