## Reviewer helper for lane wt-famlink: library paths for either arm.
## source() it with ARM set to "lane" or "base" beforehand.
if (!exists("ARM")) ARM <- "lane"
own <- if (identical(ARM, "lane")) "C:/Users/adf44/source/r/famlink-lib" else
  character()
.libPaths(c(own,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM", ARM, "frmtmb from", dirname(system.file(package = "frmtmb")), "\n")
