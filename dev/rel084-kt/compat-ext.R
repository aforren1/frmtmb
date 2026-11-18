# Merge step 6 of the 0.84.0 brief: test-compat.R on rellib-r23, with no
# extension loaded and with all eight loaded first (their registrations
# add frm_compat() rows).
#   Rscript dev/rel084-kt/compat-ext.R none|all
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
arm <- commandArgs(trailingOnly = TRUE)[1L]
suppressPackageStartupMessages(library(frmtmb))
ext <- paste0("frmtmb.", c("coupling", "eam", "latent", "learn", "ode",
                           "sample", "spatial", "survival"))
if (identical(arm, "all")) {
  for (p in ext) suppressPackageStartupMessages(library(p, character.only = TRUE))
}
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("loaded:", paste(intersect(ext, loadedNamespaces()), collapse = " "), "\n")
f <- "C:/Users/adf44/source/r/frmtmb-wt-rel084/tests/testthat/test-compat.R"
r <- as.data.frame(testthat::test_file(f, package = "frmtmb",
                                       env = testthat::test_env("frmtmb"),
                                       reporter = "silent"))
cat("RESULT", arm, "pass=", sum(r$passed), "fail=", sum(r$failed),
    "err=", sum(r$error), "skip=", sum(r$skipped), "warn=", sum(r$warning), "\n")
