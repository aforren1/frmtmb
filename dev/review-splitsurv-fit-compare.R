# Compare dev/splitsurv-fit-dump.R's two arms: Rscript dev/splitsurv-fit-compare.R
b <- readRDS("dev/review-splitsurv-out/fit-before.rds")
a <- readRDS("dev/review-splitsurv-out/fit-after.rds")
stopifnot(identical(names(a), names(b)))
for (nm in names(b)) {
  x <- a[[nm]]; y <- b[[nm]]
  if (is.list(y) && !is.data.frame(y)) {
    flds <- setdiff(names(y), "warnings")
    same <- vapply(flds, function(k) identical(x[[k]], y[[k]]), NA)
    cat(sprintf("%-18s identical: %s%s; warnings identical: %s\n", nm,
                all(same), if (!all(same)) paste0(" (differs: ",
                paste(flds[!same], collapse = ","), ")") else "",
                identical(x$warnings, y$warnings)))
  } else {
    cat(sprintf("%-18s identical: %s\n", nm, identical(x, y)))
  }
}
