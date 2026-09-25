# Library setup for every phase3a script. PHASE3A_ARM chooses the build:
# "base" reads the released 0.62.0 build in rellib-r3 only; anything else
# puts this lane's private library first.
PHASE3A_LIB <- "C:/Users/adf44/source/r/phase3a-lib"
PHASE3A_BASE <- "C:/Users/adf44/source/r/rellib-r3"
PHASE3A_USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
phase3a_arm <- Sys.getenv("PHASE3A_ARM", "lane")
if (identical(phase3a_arm, "base")) {
  .libPaths(c(PHASE3A_BASE, PHASE3A_USER))
} else {
  .libPaths(c(PHASE3A_LIB, PHASE3A_BASE, PHASE3A_USER))
}
phase3a_where <- function(pkg) {
  cat(sprintf("%s %s from %s\n", pkg, format(utils::packageVersion(pkg)),
              dirname(find.package(pkg))))
}
