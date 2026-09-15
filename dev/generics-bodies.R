# Does any generic in the ownership table carry WORK in its body?
#
# `hypothesis()` did: it armed the reserved-name shadowing note around
# `UseMethod()`.  Anything a replaceable generic does is simply not
# done once the binding resolves to the owner's generic, which is a
# bare `UseMethod()`.  That cost 8 assertions under R CMD check and
# nothing in a one-file-per-process run, so it needs a guard rather
# than a memory.
#
# Run with NO owner loaded, so every active binding hands back its own
# fallback.
LIB <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
tab <- get("frm_generic_owners", envir = ns)
bad <- character()
for (g in names(tab)) {
  f <- tryCatch(get(g), error = function(e) NULL)
  if (!is.function(f)) { bad <- c(bad, paste0(g, ": absent")); next }
  b <- body(f)
  # a single `{ UseMethod(...) }` carries no work; more than one
  # statement inside the braces does
  if (is.call(b) && identical(b[[1L]], as.name("{"))) {
    b <- if (length(b) == 2L) b[[2L]] else b
  }
  ok <- is.call(b) && identical(b[[1L]], as.name("UseMethod"))
  if (!ok) {
    bad <- c(bad, sprintf("%s: %s", g,
                          paste(trimws(deparse(body(f))),
                                collapse = " ")))
  }
}
cat(sprintf("generics in the table: %d\n", length(tab)))
cat(sprintf("carrying work in the body: %d\n", length(bad)))
for (b in bad) cat("  ", b, "\n")
cat("owners loaded during this check:",
    paste(intersect(unique(unlist(tab)), loadedNamespaces()),
          collapse = ",") , "\n")
cat("DONE\n")
