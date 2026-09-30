# Reviewer: compare the two arms of formula2-rev-regress.R element by element
a <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-formula2/dev/formula2-rev-regress-before.rds")
b <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-formula2/dev/formula2-rev-regress-after.rds")
stopifnot(identical(names(a), names(b)))
nd <- 0; ne <- 0; nerr <- 0
for (nm in names(a)) for (el in union(names(a[[nm]]), names(b[[nm]]))) {
  same <- identical(a[[nm]][[el]], b[[nm]][[el]])
  isErr <- is.list(a[[nm]][[el]]) && inherits(a[[nm]][[el]][["value"]], "ERR")
  if (isErr) nerr <- nerr + 1
  if (same) ne <- ne + 1 else {
    nd <- nd + 1
    cat("DIFF", nm, el, "\n")
    print(all.equal(a[[nm]][[el]], b[[nm]][[el]]))
  }
}
cat(sprintf("COMPARE cases=%d elements identical=%d differ=%d (errors in both arms=%d)\n",
            length(a), ne, nd, nerr))
for (nm in names(a)) for (el in names(a[[nm]]))
  if (is.list(a[[nm]][[el]]) && inherits(a[[nm]][[el]][["value"]], "ERR"))
    cat("ERRBOTH", nm, el, ":", substr(a[[nm]][[el]]$value, 1, 110), "\n")
