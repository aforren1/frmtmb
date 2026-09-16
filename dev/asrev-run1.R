## One test FILE, one R process. Reviewer's copy.
## Usage: Rscript dev/asrev-run1.R <pkgdir> <filter> [lib]
##   lib defaults to the reviewer library; pass "base" for the shared
##   reference library alone (0.57.0), or "mixed" for the reviewer core
##   with the BASE extensions behind it (the ordering-constraint case).
args <- commandArgs(trailingOnly = TRUE)
pkgdir <- args[1]
filt <- args[2]
mode <- if (length(args) >= 3L) args[3] else "own"
own <- "C:/Users/adf44/source/r/asrev-lib"
paths <- switch(mode,
  own  = c(own, "C:/Users/adf44/source/r/pinlib",
           "C:/Users/adf44/AppData/Local/R/win-library/4.6"),
  base = c("C:/Users/adf44/source/r/rellib-r3",
           "C:/Users/adf44/source/r/pinlib",
           "C:/Users/adf44/AppData/Local/R/win-library/4.6"),
  # core + frmtmb.sample from the change, every OTHER extension from the
  # 0.57.0 reference: exactly what a consolidating session sees if it
  # installs two packages and stops
  mixed = c("C:/Users/adf44/source/r/asrev-mixed",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
.libPaths(paths)
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n")
if (!identical(pkgdir, ".")) {
  pkg <- basename(pkgdir)
  suppressMessages(library(pkg, character.only = TRUE))
  cat(pkg, "from:", dirname(system.file(package = pkg)), "\n")
}
for (e in c("frmtmb.latent", "frmtmb.spline", "frmtmb.coupling")) {
  p <- system.file(package = e)
  if (nzchar(p)) cat(e, "would load from:", dirname(p), "\n")
}
this_pkg <- if (identical(pkgdir, ".")) "frmtmb" else basename(pkgdir)
setwd(file.path(pkgdir, "tests", "testthat"))
t0 <- Sys.time()
res <- testthat::test_dir(".", filter = filt, package = this_pkg,
                          reporter = "summary", stop_on_failure = FALSE)
d <- as.data.frame(res)
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("\nASREV %s pass %d fail %d error %d skip %d warn %d in %.0f s\n",
            filt, sum(d$passed), sum(d$failed), sum(d$error),
            sum(d$skipped), sum(d$warning), secs))
