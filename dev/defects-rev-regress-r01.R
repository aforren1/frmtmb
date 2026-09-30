# Reviewer of lane defects, recheck: which regression items moved between
# the lane build of the first review and the punch-round build.
a <- readRDS("dev/defects-rev-log/regress-lane-round0.rds")
b <- readRDS("dev/defects-rev-log/regress-lane.rds")
for (m in names(a)) for (k in names(a[[m]])) {
  if (!identical(a[[m]][[k]], b[[m]][[k]])) {
    cat("MOVED", m, k, " NA before:", sum(is.na(a[[m]][[k]])), " after:",
        sum(is.na(b[[m]][[k]])), "\n")
  }
}
