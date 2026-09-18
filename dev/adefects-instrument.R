# The instrument check: is the library every number was measured on
# built from THIS worktree's source?
#
#   Rscript dev/adefects-instrument.R > dev/adefects-log/instrument.txt 2>&1
#
# Every top-level function in each package's R/ is deparsed and compared
# against the installed namespace's. A stale install is the failure mode
# that makes a whole lane's measurements meaningless, and a version
# string cannot see it: nothing in this lane changes a version.
.libPaths(c("C:/Users/adf44/source/r/adefects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})

chk <- function(pkg, dir) {
  ns <- asNamespace(pkg)
  e <- new.env()
  for (f in list.files(dir, "[.]R$", full.names = TRUE)) {
    sys.source(f, envir = e)
  }
  nm <- ls(e, all.names = TRUE)
  nm <- nm[vapply(nm, function(n) is.function(get(n, envir = e)), NA)]
  same <- 0L
  diff <- character(0)
  for (n in nm) {
    if (!exists(n, envir = ns, inherits = FALSE)) {
      diff <- c(diff, paste0(n, " (absent from the namespace)"))
      next
    }
    a <- paste(deparse(get(n, envir = e)), collapse = "\n")
    b <- paste(deparse(get(n, envir = ns, inherits = FALSE)), collapse = "\n")
    if (identical(a, b)) same <- same + 1L else diff <- c(diff, n)
  }
  cat(pkg, ":", same, "of", length(nm),
      "functions identical to the source\n")
  if (length(diff)) {
    cat("  differing:", paste(utils::head(diff, 20), collapse = ", "), "\n")
  }
  # the inverse case, so that "identical" is not identical to anything:
  # a function with one line changed must be reported as differing
  invisible(length(diff) == 0L)
}

ok1 <- chk("frmtmb", "R")
ok2 <- chk("frmtmb.sample", "extensions/frmtmb.sample/R")
cat("frmtmb versions:", format(packageVersion("frmtmb")),
    format(packageVersion("frmtmb.sample")), "\n")
if (!(ok1 && ok2)) stop("the installed build is not this worktree's source")
cat("the lane library is this worktree's source\n")
