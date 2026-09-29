## Confirm the two arms of rev-13 really load different core builds.
arm <- Sys.getenv("FRMTMB_LIB", "lane")
core <- if (identical(arm, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(core, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("arm =", arm, "\n")
cat("  frmtmb path        =", getNamespaceInfo("frmtmb", "path"), "\n")
cat("  frmtmb.sample path =", getNamespaceInfo("frmtmb.sample", "path"), "\n")
cat("  has thres_pin_of_fit =",
    exists("thres_pin_of_fit", envir = asNamespace("frmtmb"),
           inherits = FALSE), "\n")
cat("  prior_entry_label has the length-1 branch =",
    any(grepl("length\\(e\\$idx\\) != 1L",
              deparse(frmtmb:::prior_entry_label))), "\n")
