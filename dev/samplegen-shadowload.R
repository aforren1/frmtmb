# The exposure the two new twin directives add. R resolves
# `S3method(owner::generic, class)` against the owner: at once when the
# owner is ALREADY loaded, from the owner's load hook when it loads
# later. So an owner on the library path that lacks the generic can
# break loading. Core's lane found this for `posterior`
# (dev/generics-findings.md, defect 3). This builds the two cases this
# lane's directives add, in both orders.
#
# Each shadow exports every name the real package exports, as a stub
# generic, EXCEPT the one name, so that nothing else frmtmb or
# frmtmb.sample registers on it is what fails.
#
#   Rscript dev/samplegen-shadowload.R <LIB> <posterior|gratia> [after]
# With `after`, frmtmb.sample is loaded FIRST and the shadow second.
av <- commandArgs(trailingOnly = TRUE)
after <- length(av) >= 3L && identical(av[3], "after")
.libPaths(c(av[1], "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkg <- av[2]
drop <- c(posterior = "rhat", gratia = "posterior_samples")[[pkg]]
ex <- parseNamespaceFile(pkg, dirname(find.package(pkg)))$exports
ex <- setdiff(ex, drop)
ex <- ex[grepl("^[A-Za-z.][A-Za-z0-9._]*$", ex)]
lib <- file.path(tempdir(), paste0("samplegen-shadow-", pkg))
d <- file.path(lib, pkg)
dir.create(file.path(d, "R"), recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste0("Package: ", pkg), "Version: 99.0", "Title: t",
             "Author: t", "Maintainer: t <t@t.tt>", "Description: t.",
             "License: GPL-2", "Built: R 4.6.1; ; ; windows"),
           file.path(d, "DESCRIPTION"))
writeLines(sprintf("export(%s)", ex), file.path(d, "NAMESPACE"))
writeLines(sprintf("`%s` <- function(x, ...) UseMethod(\"%s\")", ex, ex),
           file.path(d, "R", pkg))
.libPaths(c(lib, .libPaths()))
try_load <- function(p) {
  tryCatch({ suppressMessages(loadNamespace(p)); "loaded" },
           error = function(e) paste("FAILED:", conditionMessage(e)))
}
cat("BUILD", basename(av[1]), if (after) "sample first" else "shadow first",
    "\n")
if (after) {
  cat("  FRMTMB.SAMPLE", substr(try_load("frmtmb.sample"), 1, 100), "\n")
  cat("  SHADOW", pkg, "without", drop, substr(try_load(pkg), 1, 100), "\n")
} else {
  cat("  SHADOW", pkg, "without", drop, try_load(pkg), "\n")
  cat("  FRMTMB.SAMPLE", substr(try_load("frmtmb.sample"), 1, 100), "\n")
}
