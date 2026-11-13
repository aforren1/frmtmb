# Lane gpstop: shared helpers. Arm "lane" reads wt-gpstop-lib over
# rellib-r20; arm "base" reads rellib-r20 alone (frmtmb 0.81.0 at e8abed0a).
# The names of lanes fitedge's and arstep's helpers are kept, so their
# scripts run here with the source line changed.
GS_WT <- "C:/Users/adf44/source/r/frmtmb-wt-gpstop"
GS_LIB <- "C:/Users/adf44/source/r/wt-gpstop-lib"
GS_BASE <- "C:/Users/adf44/source/r/rellib-r20"
GS_USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
GS_OUT <- file.path(GS_WT, "dev/gpstop-out")
dir.create(GS_OUT, showWarnings = FALSE)
FE_WT <- GS_WT
RV_OUT <- GS_OUT
REV_OUT <- GS_OUT
gs_libs <- function(arm = c("lane", "base")) {
  arm <- match.arg(arm)
  .libPaths(c(if (arm == "lane") GS_LIB, GS_BASE, GS_USER))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
             FRMTMB_STAN_CACHE = file.path(GS_WT, "dev/stan-cache"),
             NOT_CRAN = "true")
  suppressPackageStartupMessages(library(frmtmb))
  cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
      as.character(packageVersion("frmtmb")), "arm", arm, "\n")
  invisible(arm)
}
fe_libs <- function(arm = "lane") gs_libs(if (arm == "base") "base" else "lane")
rev_libs <- fe_libs
rv_libs <- fe_libs

# A fit with its warnings kept (class and text), never an error thrown.
fe_fit <- function(...) {
  w <- character(0)
  f <- tryCatch(withCallingHandlers(frm(...), warning = function(e) {
    w <<- c(w, paste0(class(e)[1], ": ", conditionMessage(e)))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(fit = f, warn = w)
}
rev_fit <- function(fml, data, ...) {
  w <- character(0)
  fit <- withCallingHandlers(
    tryCatch(frmtmb::frm(fml, data = data, ...), error = function(e) e),
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, warn = w)
}
rv_put <- function(name, f) {
  ns <- asNamespace("frmtmb")
  environment(f) <- ns
  unlockBinding(name, ns); assign(name, f, envir = ns); lockBinding(name, ns)
}
rv_quiet <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  list(value = v, warnings = w)
}
