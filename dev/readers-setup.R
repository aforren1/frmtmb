# Library setup for the readers lane's scripts. READERS_ARM=before reads the
# base build rellib-r18 alone; anything else puts the lane library first.
arm <- Sys.getenv("READERS_ARM", "after")
lib_lane <- "C:/Users/adf44/source/r/wt-readers-lib"
lib_base <- "C:/Users/adf44/source/r/rellib-r18"
lib_user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "before")) c(lib_base, lib_user) else
  c(lib_lane, lib_base, lib_user))
suppressPackageStartupMessages(library(frmtmb))
cat("arm:", arm, " lib:", dirname(find.package("frmtmb")),
    " frmtmb", format(packageVersion("frmtmb")), "\n")
