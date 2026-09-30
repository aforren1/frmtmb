# Reviewer (re-check after punch round 1), lane ordinal: does the lane library hold the worktree source?
# Installs a copy of the worktree (core and frmtmb.sample) into a scratch
# library under the session scratchpad, then compares every namespace
# object of that install with the lane library's, by deparse.
# Output: dev/ordinal-rev2-log-libmatch.txt
sp <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad"
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ordinal"
scr <- file.path(sp, "ordrev-lib2")
lane <- "C:/Users/adf44/source/r/wt-ordinal-lib"
src <- file.path(sp, "ordrev-src2")
unlink(src, recursive = TRUE)
dir.create(src)
for (p in c("core", "sample")) {
  from <- if (p == "core") wt else file.path(wt, "extensions/frmtmb.sample")
  to <- file.path(src, p)
  dir.create(to)
  keep <- setdiff(list.files(from, all.files = FALSE, no.. = TRUE),
                  c("dev", "extensions", "docs", "vignettes"))
  file.copy(file.path(from, keep), to, recursive = TRUE)
}
R <- file.path(R.home("bin"), "R.exe")
for (p in c("core", "sample")) {
  st <- system2(R, c("CMD", "INSTALL", paste0("--library=", scr),
                     "--no-multiarch", "--no-test-load", file.path(src, p)),
                stdout = TRUE, stderr = TRUE,
                env = paste0("R_LIBS=", paste(c(scr, lane,
                  "C:/Users/adf44/source/r/rellib-r4",
                  "C:/Users/adf44/AppData/Local/R/win-library/4.6"),
                  collapse = ";")))
  cat(tail(st, 2), sep = "\n")
}
cmp <- function(pkg) {
  a <- new.env(); b <- new.env()
  lazyLoad(file.path(scr, pkg, "R", pkg), envir = a)
  lazyLoad(file.path(lane, pkg, "R", pkg), envir = b)
  na <- sort(ls(a, all.names = TRUE)); nb <- sort(ls(b, all.names = TRUE))
  cat(pkg, ": objects scratch", length(na), "lane", length(nb),
      " only-scratch:", paste(setdiff(na, nb), collapse = ","),
      " only-lane:", paste(setdiff(nb, na), collapse = ","), "\n")
  diff <- character(0)
  for (n in intersect(na, nb)) {
    x <- get(n, a); y <- get(n, b)
    dx <- if (is.function(x)) deparse(x) else deparse(x)
    dy <- if (is.function(y)) deparse(y) else deparse(y)
    if (!identical(dx, dy)) diff <- c(diff, n)
  }
  cat(pkg, ": objects that differ:", length(diff),
      if (length(diff)) paste(head(diff, 30), collapse = ","), "\n")
  ns_a <- readLines(file.path(scr, pkg, "NAMESPACE"))
  ns_b <- readLines(file.path(lane, pkg, "NAMESPACE"))
  cat(pkg, ": NAMESPACE identical:", identical(ns_a, ns_b), "\n")
}
sink(file.path(wt, "dev/ordinal-rev2-log-libmatch.txt"), split = TRUE)
cmp("frmtmb")
cmp("frmtmb.sample")
sink()
