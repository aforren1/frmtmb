# Reviewer re-check: identical() of the downstream readers, base vs lane.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
b <- readRDS(file.path(wt, "dev/ordinal-rev2-disc-base.rds"))
l <- readRDS(file.path(wt, "dev/ordinal-rev2-disc-lane.rds"))
for (k in names(b)) cat(sprintf("%-12s identical: %s\n", k, identical(b[[k]], l[[k]])))
