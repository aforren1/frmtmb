# Reviewer 2 prelude. LANE = phase3b-lib first; BASE = rellib-r3.
r2_lib <- function(which = c("lane", "base")) {
  which <- match.arg(which)
  first <- if (which == "lane") "C:/Users/adf44/source/r/phase3b-lib" else
    "C:/Users/adf44/source/r/rellib-r3"
  .libPaths(c(first, "C:/Users/adf44/source/r/pinlib",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  invisible(.libPaths())
}
