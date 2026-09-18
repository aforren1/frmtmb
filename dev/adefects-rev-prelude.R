# Reviewer prelude. ARM is "lane" or "base".
USERLIB <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
LANELIB <- "C:/Users/adf44/source/r/adefects-lib"
BASELIB <- "C:/Users/adf44/source/r/rellib-r3"
arm <- Sys.getenv("REV_ARM", "lane")
if (identical(arm, "lane")) {
  .libPaths(c(LANELIB, BASELIB, USERLIB))
} else {
  .libPaths(c(BASELIB, USERLIB))
}
WT <- "C:/Users/adf44/source/r/frmtmb-wt-adefects"
cat("ARM", arm, "\n")
cat("frmtmb", as.character(utils::packageVersion("frmtmb")),
    "at", dirname(system.file(package = "frmtmb")), "\n")
