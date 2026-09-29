## Reviewer: confirm the library each arm actually loads, and that the
## worker's installed build carries the worktree's new code verbatim.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib   =", lib, "\n")
cat("found =", dirname(dirname(getNamespaceInfo("frmtmb", "path"))), "\n")
cat("ver   =", as.character(utils::packageVersion("frmtmb")), "\n")
has <- exists("thres_pin_of_fit", envir = asNamespace("frmtmb"),
              inherits = FALSE)
cat("thres_pin_of_fit present =", has, "\n")
if (has) {
  ## the installed body against the worktree source, so a stale install
  ## cannot pass as the change under review
  src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-thresrefit/R/thres.R")
  i0 <- grep("^thres_pin_of_fit <- function", src)
  i1 <- grep("^thres_pin_apply <- function", src)
  cat("source lines: of_fit at", i0, " apply at", i1, "\n")
  f <- get("thres_pin_of_fit", envir = asNamespace("frmtmb"))
  g <- get("thres_pin_apply", envir = asNamespace("frmtmb"))
  cat("--- installed thres_pin_of_fit ---\n")
  cat(paste(deparse(f), collapse = "\n"), "\n")
  cat("--- installed thres_pin_apply ---\n")
  cat(paste(deparse(g), collapse = "\n"), "\n")
  cat("assemble_frame formals:",
      paste(names(formals(get("assemble_frame",
                              envir = asNamespace("frmtmb")))),
            collapse = ", "), "\n")
  cat("influence body has thres_pin:",
      any(grepl("thres_pin", deparse(frmtmb:::influence.frmtmb_fit))), "\n")
  cat("draw_prior_entry has the refusal:",
      any(grepl("parameters at once",
                deparse(frmtmb:::draw_prior_entry))), "\n")
  cat("prior_entry_label has the length-1 branch:",
      any(grepl("length\\(e\\$idx\\) != 1L",
                deparse(frmtmb:::prior_entry_label))), "\n")
}
cat("OK\n")
