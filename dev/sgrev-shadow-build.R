# Build the partial-owner shadows in a process of their own. The first
# spelling read the real package's exports IN the loading process, so
# the real owner was already loaded when the shadow was asked for and
# the guard "the shadow really lacks the name" read FALSE. Building and
# loading must not share a process.
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
SHADOWLIB <- "C:/Users/adf44/source/r/sgrev-shadowlib"
unlink(SHADOWLIB, recursive = TRUE)
dir.create(SHADOWLIB, recursive = TRUE, showWarnings = FALSE)
for (owner in c("posterior", "gratia")) {
  drop <- if (identical(owner, "posterior")) "rhat" else
    "posterior_samples"
  # exports read from the installed NAMESPACE file, which needs no
  # load, so the real package is never in this process's way
  lib <- dirname(system.file(package = owner))
  info <- parseNamespaceFile(owner, lib)
  real <- unique(c(info$exports,
                   unlist(lapply(info$exportPatterns, function(p) NULL))))
  keep <- setdiff(real, drop)
  # Non-syntactic exports such as `%in%` and `draws_of<-` are dropped:
  # R's namespaceExport() reports them undefined in a package this
  # small, which is a quirk of the shadow and not of the thing under
  # test. Measured on a four-name package in dev/sgrev-out/.
  keep <- keep[grepl("^[.a-zA-Z][._a-zA-Z0-9]*$", keep)]
  stopifnot(drop %in% real, !(drop %in% keep), length(keep) > 20)
  src <- file.path(SHADOWLIB, paste0("src-", owner))
  dir.create(file.path(src, "R"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(paste0("Package: ", owner), "Version: 0.0.1",
               "Title: shadow", "Description: A shadow.",
               "Author: x", "Maintainer: x <x@x.com>", "License: GPL-2"),
             file.path(src, "DESCRIPTION"))
  writeLines(paste0("export(", paste0("`", keep, "`", collapse = ","),
                    ")"), file.path(src, "NAMESPACE"))
  writeLines(paste0(vapply(keep, function(n)
    paste0("`", n, "` <- function(...) NULL"), ""), collapse = "\n"),
    file.path(src, "R", "a.R"))
  o <- system2(file.path(R.home("bin"), "R"),
               c("CMD", "INSTALL", paste0("--library=", shQuote(SHADOWLIB)),
                 "--no-multiarch", "--no-docs", "--no-byte-compile",
                 "--no-test-load", shQuote(src)),
               stdout = TRUE, stderr = TRUE)
  cat("shadow ", owner, " without ", drop, ": ",
      if (any(grepl("DONE", o))) "installed" else "FAILED",
      "  (", length(keep), " exports kept)\n", sep = "")
  unlink(src, recursive = TRUE)
}
cat("DONE\n")
