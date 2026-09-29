# Reviewer: identical() on the emmeans output of fits with no
# group-indexed smooth, between the two arms.
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
a <- readRDS("dev/resmooth-rev-emm-lane.rds")
b <- readRDS("dev/resmooth-rev-emm-base.rds")
for (nm in names(a)) {
  id <- identical(a[[nm]], b[[nm]])
  md <- NA_real_
  if (is.data.frame(a[[nm]]) && is.data.frame(b[[nm]])) {
    na <- vapply(a[[nm]], is.numeric, NA)
    md <- max(abs(unlist(a[[nm]][na]) - unlist(b[[nm]][na])))
  }
  cat(sprintf("%-14s identical=%-5s max|numeric diff|=%.3e\n", nm, id, md))
}
