od <- "C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-rev-out"
b <- readRDS(file.path(od, "ord2-base.rds")); l <- readRDS(file.path(od, "ord2-lane.rds"))
for (nm in c("influence_unused", "influence_ofactor")) {
  x <- b[[nm]]$value; y <- l[[nm]]$value
  cat(nm, " warnings identical:", identical(b[[nm]]$warnings, l[[nm]]$warnings), "\n")
  for (k in names(x)) if (!identical(x[[k]], y[[k]])) {
    cat("  element", k, "differs; class", class(x[[k]])[1], "\n")
    if (is.numeric(x[[k]])) cat("    max abs diff", max(abs(x[[k]] - y[[k]])), " max rel",
      max(abs(x[[k]] - y[[k]]) / pmax(abs(x[[k]]), 1e-300)), "\n") else str(x[[k]], max.level = 1)
  }
}
