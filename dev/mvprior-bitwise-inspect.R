# Lane wt-mvprior: what differs in a bitwise entry, printed side by side.
# Usage: Rscript dev/mvprior-bitwise-inspect.R "<design key>" <quantity>
a <- commandArgs(trailingOnly = TRUE)
OUT <- "C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-log"
x <- readRDS(file.path(OUT, "bitwise-base.rds"))[[a[1]]][[a[2]]]
y <- readRDS(file.path(OUT, "bitwise-lane.rds"))[[a[1]]][[a[2]]]
cat("== base\n"); str(x, max.level = 3)
cat("== lane\n"); str(y, max.level = 3)
print(all.equal(x, y))
# for a list of prior entries: the same set in another order?
if (identical(a[2], "entries") && is.list(x) && is.list(y)) {
  key <- function(e) paste(e$comp, paste(e$idx, collapse = ","))
  ox <- x[order(vapply(x, key, ""))]
  oy <- y[order(vapply(y, key, ""))]
  cat("same entries, order aside:", identical(ox, oy), "\n")
}
