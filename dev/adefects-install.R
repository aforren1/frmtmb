# Roxygenise and install this worktree's packages into the lane's
# PRIVATE library only. Every other library on this machine is read
# only to this lane, including the round's reference build rellib-r3,
# which holds the base commit f8b45ef and is what the numbers are
# measured against.
#
#   Rscript dev/adefects-install.R [frmtmb] [frmtmb.sample] ...
lib <- "C:/Users/adf44/source/r/adefects-lib"
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(identical(normalizePath(.libPaths()[1], winslash = "/"),
                    normalizePath(lib, winslash = "/")))
root <- "C:/Users/adf44/source/r/frmtmb-wt-adefects"
pkgs <- commandArgs(trailingOnly = TRUE)
if (!length(pkgs)) pkgs <- c("frmtmb", "frmtmb.sample")
# Roxygenise in a CHILD process, one per package. roxygen2 loads the
# package with pkgload BEFORE it rewrites NAMESPACE, so a namespace
# loaded earlier in this session is the one an extension's
# `importFrom(frmtmb, ...)` resolves against: with both packages in one
# process, frmtmb.sample stopped at "'log_lik' is not an exported object
# from 'namespace:frmtmb'" on a frmtmb that exports it.
lp <- paste0(deparse(.libPaths()), collapse = "")
for (p in pkgs) {
  path <- if (p == "frmtmb") root else file.path(root, "extensions", p)
  f <- tempfile(fileext = ".R")
  writeLines(c(sprintf(".libPaths(%s)", lp),
               sprintf("roxygen2::roxygenise(%s)", shQuote(path))), f)
  st <- system2(file.path(R.home("bin"), "Rscript.exe"),
                c("--vanilla", shQuote(f)))
  unlink(f)
  if (st != 0) stop("roxygenise of ", p, " failed")
  st <- system2(file.path(R.home("bin"), "R.exe"),
                c("CMD", "INSTALL", "--no-multiarch",
                  paste0("--library=", shQuote(lib)), shQuote(path)))
  if (st != 0) stop("install of ", p, " failed")
}
cat("installed", pkgs, "into", lib, "\n")
