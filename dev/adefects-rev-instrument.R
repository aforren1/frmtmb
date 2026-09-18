source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")

check_pkg <- function(pkg, srcdir) {
  ns <- asNamespace(pkg)
  files <- list.files(srcdir, pattern = "[.]R$", full.names = TRUE)
  n_ok <- 0L; n_bad <- 0L; bad <- character()
  n_missing <- 0L; miss <- character()
  for (f in files) {
    ex <- parse(f, keep.source = FALSE)
    for (e in ex) {
      if (!(is.call(e) && length(e) == 3L &&
            as.character(e[[1L]]) %in% c("<-", "=") &&
            is.name(e[[2L]]) && is.call(e[[3L]]) &&
            identical(e[[3L]][[1L]], as.name("function")))) next
      nm <- as.character(e[[2L]])
      if (!exists(nm, envir = ns, inherits = FALSE)) {
        n_missing <- n_missing + 1L; miss <- c(miss, nm); next
      }
      if (bindingIsActive(nm, ns)) next
      inst <- get(nm, envir = ns)
      if (!is.function(inst)) next
      src <- eval(e[[3L]], envir = globalenv())
      a <- paste(deparse(src), collapse = "\n")
      b <- paste(deparse(inst), collapse = "\n")
      if (identical(a, b)) n_ok <- n_ok + 1L else {
        n_bad <- n_bad + 1L; bad <- c(bad, nm)
      }
    }
  }
  cat(sprintf("%s: ok=%d mismatch=%d missing=%d\n", pkg, n_ok, n_bad,
              n_missing))
  if (n_bad) cat("  MISMATCH:", paste(head(bad, 20), collapse = ", "), "\n")
  if (n_missing) cat("  MISSING:", paste(head(miss, 20), collapse = ", "), "\n")
}

suppressMessages(library(frmtmb))
check_pkg("frmtmb", file.path(WT, "R"))
suppressMessages(library(frmtmb.sample))
cat("frmtmb.sample", as.character(utils::packageVersion("frmtmb.sample")),
    "\n")
check_pkg("frmtmb.sample", file.path(WT, "extensions/frmtmb.sample/R"))
