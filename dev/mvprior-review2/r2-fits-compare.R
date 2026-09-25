root <- "C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2"
a <- readRDS(file.path(root, "r2-fits-lane.rds")); b <- readRDS(file.path(root, "r2-fits-base.rds"))
tot <- 0; same <- 0
for (nm in intersect(names(a), names(b))) {
  if (!is.list(a[[nm]]) || !is.list(b[[nm]])) { cat(nm, ": lane", class(a[[nm]])[1], "base", class(b[[nm]])[1], "\n"); next }
  for (q in intersect(names(a[[nm]]), names(b[[nm]]))) {
    tot <- tot + 1; s <- identical(a[[nm]][[q]], b[[nm]][[q]]); same <- same + s
    if (!s) { cat(nm, q, "DIFF\n"); if (is.data.frame(a[[nm]][[q]])) { print(a[[nm]][[q]]); print(b[[nm]][[q]]) } }
  }
}
cat("quantities", tot, "identical", same, "\n")
