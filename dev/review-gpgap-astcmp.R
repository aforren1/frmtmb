# Short check of fix round 1, item 6: does the lane's installed build hold the code
# of the worktree as it stands? Every object of each namespace, srcrefs dropped,
# against a fresh install of the current sources. And the lane's last edit scripts
# (ed41.R, ed42.R): is every changed line a comment or roxygen?
#   Rscript dev/review-gpgap-astcmp.R <pkg> <lib a> <lib b>
a <- commandArgs(trailingOnly = TRUE)
pkg <- a[1]
read_ns <- function(lib) {
  e <- new.env()
  lazyLoad(file.path(lib, pkg, "R", pkg), envir = e)
  e
}
A <- read_ns(a[2]); B <- read_ns(a[3])
strip <- function(x) {
  if (is.function(x)) {
    x <- utils::removeSource(x)
    return(list(formals(x), body(x)))
  }
  if (is.environment(x)) return("env")
  x
}
na <- sort(ls(A, all.names = TRUE)); nb <- sort(ls(B, all.names = TRUE))
cat(sprintf("%s: %d objects in %s, %d in %s; only one side: %s\n", pkg, length(na),
            basename(a[2]), length(nb), basename(a[3]),
            paste(c(setdiff(na, nb), setdiff(nb, na)), collapse = ",")))
diffs <- character()
for (n in intersect(na, nb)) {
  if (!isTRUE(all.equal(strip(A[[n]]), strip(B[[n]]), check.environment = FALSE))) {
    if (!identical(deparse(strip(A[[n]])), deparse(strip(B[[n]])))) diffs <- c(diffs, n)
  }
}
cat("objects whose code differs:", length(diffs), paste(utils::head(diffs, 20), collapse = ","), "\n")
