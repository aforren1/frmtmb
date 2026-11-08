# Shared by lane fitedge's scripts. Arm "lane" reads the lane's private
# library first; arm "base" reads rellib-r18 (frmtmb 0.79.1 at 3cd2a0c8).
FE_WT <- "C:/Users/adf44/source/r/frmtmb-wt-fitedge"
FE_LIB <- "C:/Users/adf44/source/r/wt-fitedge-lib"
FE_BASE <- "C:/Users/adf44/source/r/rellib-r18"
FE_USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
FE_OUT <- file.path(FE_WT, "dev/fitedge-out")
dir.create(FE_OUT, showWarnings = FALSE)
fe_libs <- function(arm = c("lane", "base", "dev")) {
  arm <- match.arg(arm)
  .libPaths(c(if (arm == "lane") FE_LIB, FE_BASE, FE_USER))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
             FRMTMB_STAN_CACHE = file.path(FE_WT, "dev/stan-cache"),
             NOT_CRAN = "true")
  # "dev" reads the worktree's R/ unbuilt, for a probe while the lane's
  # library is in use by a running process; no claim rests on it
  if (arm == "dev") {
    suppressPackageStartupMessages(pkgload::load_all(FE_WT, quiet = TRUE,
                                                     export_all = FALSE))
    cat("lib: load_all", FE_WT, "\n")
    return(invisible(arm))
  }
  suppressPackageStartupMessages(library(frmtmb))
  cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
      as.character(packageVersion("frmtmb")), "arm", arm, "\n")
  invisible(arm)
}

# A fit with its warnings kept, never an error thrown.
fe_fit <- function(...) {
  w <- character(0)
  f <- tryCatch(withCallingHandlers(frm(...), warning = function(e) {
    w <<- c(w, paste0(class(e)[1], ": ", conditionMessage(e)))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(fit = f, warn = w)
}
