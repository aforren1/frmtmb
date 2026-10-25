# The consolidation's new tests in test-extension-home-mixed.R, seen to
# fail without resolution 1 of dev/rel075-resolve.R: in one process on
# rellib-r13, seam_missing_parts() is put back to the structures alone
# (the merge before the resolution: smallfix's collector without
# sarfcor's residual terms), and stop_extension_parts() to a copy whose
# "session" sentence names frmtmb_register_covstruct() for every part.
# The file's in-process blocks then fail where they name the residual
# term. The child blocks read the installed build and are not affected.
#   Rscript dev/rel075-seenfail.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")),
    format(packageVersion("frmtmb")), "\n")
smp <- function(x) ns$structure_missing_parts(x)
environment(smp) <- ns
unlockBinding("seam_missing_parts", ns)
assign("seam_missing_parts", smp, envir = ns)
sep <- ns$stop_extension_parts
b <- deparse(body(sep))
b <- sub('p[["register"]] %||% "frmtmb_register_covstruct()"',
         '"frmtmb_register_covstruct()"', b, fixed = TRUE)
body(sep) <- parse(text = paste(b, collapse = "\n"))[[1]]
unlockBinding("stop_extension_parts", ns)
assign("stop_extension_parts", sep, envir = ns)
root <- "C:/Users/adf44/source/r/frmtmb-wt-rel075/tests/testthat"
env <- test_env("frmtmb")
for (h in list.files(root, "^helper-.*[.]R$", full.names = TRUE)) {
  sys.source(h, envir = env)
}
r <- as.data.frame(test_file(file.path(root, "test-extension-home-mixed.R"),
                             package = "frmtmb", env = env,
                             reporter = "silent"))
print(r[, c("test", "nb", "failed", "skipped", "error")], right = FALSE)
cat("failed in the consolidation's blocks:",
    sum(r$failed[grepl("three absent|car\\(\\) and sar\\(\\)", r$test)]),
    "\n")
