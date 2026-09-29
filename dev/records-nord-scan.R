# Scan every R source in the monorepo for a `#' @noRd` line followed by
# more roxygen prose. pkgcheck refuses that shape ("@noRd must not be
# followed by any text") and nothing in the test suite catches it, so
# this is the instrument for the rule that dev/lane-rules.md states.
#
# Usage: Rscript dev/records-nord-scan.R [root]
# Exit status 1 when it finds one, so a driver can gate on it.

a <- commandArgs(trailingOnly = TRUE)
root <- if (length(a)) a[1] else "."

dirs <- c(file.path(root, "R"),
          Sys.glob(file.path(root, "extensions", "*", "R")))
dirs <- dirs[dir.exists(dirs)]
files <- unlist(lapply(dirs, list.files, pattern = "[.]R$",
                       full.names = TRUE))

is_rox <- function(x) grepl("^[[:space:]]*#'", x)
is_nord <- function(x) grepl("^[[:space:]]*#'[[:space:]]*@noRd[[:space:]]*$", x)
# a tag line is allowed after @noRd; prose is not
is_prose <- function(x) {
  b <- sub("^[[:space:]]*#'[[:space:]]*", "", x)
  nzchar(b) && !grepl("^@", b)
}

hits <- 0L
for (f in files) {
  x <- readLines(f, warn = FALSE)
  inblock <- FALSE
  for (i in seq_along(x)) {
    if (is_nord(x[i])) { inblock <- TRUE; next }
    if (!inblock) next
    if (!is_rox(x[i])) { inblock <- FALSE; next }
    if (is_prose(x[i])) {
      hits <- hits + 1L
      cat(f, ":", i, ": ", x[i], "\n", sep = "")
    }
  }
}

cat("files scanned:", length(files), " prose-after-@noRd lines:", hits, "\n")
if (hits > 0L) quit(status = 1L)
