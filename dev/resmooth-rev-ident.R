# Reviewer: identical() between the two arms' saved fitted() output on
# fits the change must not touch, plus the gp() cost scaling.
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cmp <- function(tag) {
  a <- readRDS(file.path("dev", paste0("resmooth-rev-", tag, "-lane.rds")))
  b <- readRDS(file.path("dev", paste0("resmooth-rev-", tag, "-base.rds")))
  if (!is.list(a)) { a <- list(x = a); b <- list(x = b) }
  for (nm in names(a)) {
    id <- identical(a[[nm]], b[[nm]])
    md <- max(abs(as.vector(a[[nm]]) - as.vector(b[[nm]])))
    cat(sprintf("%-14s %-14s identical=%-5s max|diff|=%.3e\n", tag, nm,
                id, md))
  }
}
cmp("nosmooth")
cmp("nearlin")
cmp("dpar")
