# Load half of item 6: the shadow library is FIRST, so the shadow is
# the owner this process can see at all.
#   Rscript dev/sgrev-shadow-load.R <LIB> <owner> <order>
av <- commandArgs(trailingOnly = TRUE)
LIB <- av[1]; owner <- av[2]; order <- av[3]
.libPaths(c("C:/Users/adf44/source/r/sgrev-shadowlib", LIB,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
drop <- if (identical(owner, "posterior")) "rhat" else "posterior_samples"
cat("LIB ", basename(LIB), "  owner ", owner, " (no ", drop,
    ")  order ", order, "\n", sep = "")
loadsh <- function() {
  r <- tryCatch({ suppressMessages(loadNamespace(owner)); "loaded" },
                error = function(e) paste("FAILED:", conditionMessage(e)))
  cat("  SHADOW ", owner, ": ", r, "\n", sep = "")
  # the positive condition, asserted rather than assumed
  cat("  the shadow is what got loaded (it lacks ", drop, "): ",
      !(drop %in% getNamespaceExports(owner)), "\n", sep = "")
}
loadfs <- function() {
  r <- tryCatch({ suppressMessages(loadNamespace("frmtmb.sample"))
    "loaded" },
    error = function(e) paste("FAILED:", conditionMessage(e)))
  cat("  frmtmb.sample: ", r, "\n", sep = "")
}
if (identical(order, "shadow-first")) { loadsh(); loadfs() } else
  { loadfs(); loadsh() }
cat("DONE\n")
