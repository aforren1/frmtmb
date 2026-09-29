## Reviewer: identical() over two saved lists, element by element.
## usage: Rscript gradcheck-rev-15-cmp2.R <lane.rds> <base.rds>
a <- commandArgs(TRUE)
A <- readRDS(a[1]); B <- readRDS(a[2])
cat("names identical:", identical(names(A), names(B)), "\n")
ok <- TRUE
for (nm in union(names(A), names(B))) {
  x <- A[[nm]]; y <- B[[nm]]
  same <- identical(x, y)
  extra <- ""
  if (!same && is.numeric(x) && is.numeric(y) &&
      identical(dim(x), dim(y)) && length(x) == length(y)) {
    extra <- sprintf("  max|diff| = %s", format(max(abs(x - y)), digits = 6))
  }
  if (!same) ok <- FALSE
  cat(sprintf("%-16s identical: %-5s%s\n", nm, same, extra))
  if (!same && is.list(x) && is.list(y)) {
    for (k in union(names(x), names(y))) {
      if (!identical(x[[k]], y[[k]])) {
        cat("    field", k, "differs: lane=",
            paste(format(x[[k]]), collapse = ","), " base=",
            paste(format(y[[k]]), collapse = ","), "\n")
      }
    }
  }
}
cat("\nALL IDENTICAL:", ok, "\n")
