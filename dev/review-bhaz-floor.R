# Reviewer: the lane's frmtmb.survival in front of rellib-r12's frmtmb.
.libPaths(c("C:/Users/adf44/source/r/frmtmb-wt-bhaz/dev/review-bhaz-floorlib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("frmtmb", format(packageVersion("frmtmb")), "at", find.package("frmtmb"),
    "\nsurvival at", find.package("frmtmb.survival"), "\n")
r <- tryCatch(library(frmtmb.survival), error = function(e) conditionMessage(e))
cat("library(frmtmb.survival):", if (is.character(r)) r else "loaded", "\n")
