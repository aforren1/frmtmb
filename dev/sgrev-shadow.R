# Item 6. The new load failure, reproduced independently.
#   half 1: a posterior without rhat, loaded BEFORE frmtmb.sample
#   half 2: a gratia without posterior_samples, same
# Built here rather than rerun from the lane's script, and extended
# with the two questions the lane did not answer: is a partial owner
# loaded AFTER harmless, and what is the price of dropping gratia.
#
#   Rscript dev/sgrev-shadow.R <LIB> <owner> <order>
av <- commandArgs(trailingOnly = TRUE)
LIB <- av[1]; owner <- av[2]; order <- av[3]
SHADOWLIB <- file.path(tempdir(), "sgrev-shadowlib")
dir.create(SHADOWLIB, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

drop <- if (identical(owner, "posterior")) "rhat" else "posterior_samples"
# Build a package with the owner's NAME that exports everything the
# real one does EXCEPT the one name, so "a partial or older owner" is
# constructed rather than imagined.
real <- getNamespaceExports(owner)
keep <- setdiff(real, drop)
src <- file.path(tempdir(), paste0("sgrev-", owner))
unlink(src, recursive = TRUE)
dir.create(file.path(src, "R"), recursive = TRUE)
writeLines(c(paste0("Package: ", owner),
             "Version: 0.0.1", "Title: shadow", "Description: shadow.",
             "Author: x", "Maintainer: x <x@x.com>", "License: GPL-2"),
           file.path(src, "DESCRIPTION"))
writeLines(paste0("export(", paste0("`", keep, "`", collapse = ","), ")"),
           file.path(src, "NAMESPACE"))
writeLines(c(paste0(vapply(keep, function(n)
  paste0("`", n, "` <- function(...) NULL"), ""), collapse = "\n")),
  file.path(src, "R", "a.R"))
o <- system2(file.path(R.home("bin"), "R"),
             c("CMD", "INSTALL", paste0("--library=", shQuote(SHADOWLIB)),
               "--no-multiarch", "--no-docs", "--no-byte-compile",
               "--no-test-load", shQuote(src)),
             stdout = TRUE, stderr = TRUE)
cat("shadow", owner, "without", drop, ": installed",
    any(grepl("DONE", o)), "\n")
cat("  it exports", length(keep), "names and not", drop, "\n")

.libPaths(c(SHADOWLIB, .libPaths()))
loadsh <- function() {
  r <- tryCatch({ suppressMessages(loadNamespace(owner)); "loaded" },
                error = function(e) paste("FAILED:", conditionMessage(e)))
  cat("  SHADOW ", owner, ": ", r, "\n", sep = "")
  cat("  shadow really lacks ", drop, ": ",
      !(drop %in% getNamespaceExports(owner)), "\n", sep = "")
}
loadfs <- function() {
  r <- tryCatch({ suppressMessages(loadNamespace("frmtmb.sample"))
    "loaded" },
    error = function(e) paste("FAILED:", conditionMessage(e)))
  cat("  frmtmb.sample: ", r, "\n", sep = "")
}
cat("LIB ", LIB, " owner ", owner, " order ", order, "\n", sep = "")
if (identical(order, "shadow-first")) { loadsh(); loadfs() } else
  { loadfs(); loadsh() }
cat("DONE\n")
