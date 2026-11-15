# Short check of lane priors2's fix round, item 5: famend's
# test-prior-table-brms.R with the merge note's at0 block, run against a
# scratch tree (lane priors2's R/ with R/priors.R and R/fit.R three-way
# merged with famend's, git merge-file after CRLF -> LF), loaded with
# pkgload::load_all. Nothing is installed.
# Rscript priors2-short-merge.R <scratch package tree>
tree <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(
  pkgload::load_all(tree, export_all = FALSE, helpers = FALSE, quiet = TRUE))
cat("frmtmb from:", getNamespaceInfo("frmtmb", "path"), "\n")
cat("prior_start_inside:", exists("prior_start_inside",
                                  envir = asNamespace("frmtmb")),
    " open_end_rescue:", exists("open_end_rescue",
                                envir = asNamespace("frmtmb")), "\n")
f <- file.path(tree, "tests", "testthat", "test-prior-table-brms.R")
r <- as.data.frame(testthat::test_file(f, package = "frmtmb",
                                       reporter = "silent",
                                       load_package = "none"))
cat("RESULT", basename(f), "pass=", sum(r$passed), "fail=", sum(r$failed),
    "err=", sum(r$error), "skip=", sum(r$skipped), "warn=",
    sum(r$warning), "\n")
bad <- r[r$failed > 0 | r$error, c("test", "failed", "error")]
if (nrow(bad)) print(bad)
