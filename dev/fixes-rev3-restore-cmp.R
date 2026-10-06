# Reviewer of lane fixes, final check: compare dev/fixes-rev3-restore.R's
# three arms element by element, bitwise.
A <- readRDS("dev/fixes-rev3-log/restore-check.rds")
N <- readRDS("dev/fixes-rev3-log/restore-noop.rds")
B <- readRDS("dev/fixes-rev3-log/restore-base.rds")
for (k in names(A)) {
  for (el in names(A[[k]])) {
    cat(sprintf("%-8s %-9s check==noop %-5s check==base %-5s\n", k, el,
                identical(A[[k]][[el]], N[[k]][[el]]),
                identical(A[[k]][[el]], B[[k]][[el]])))
  }
}
cat("boot is a matrix:", is.matrix(A[[1]]$boot), "\n")
