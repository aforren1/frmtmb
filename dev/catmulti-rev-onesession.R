# Reviewer copy of dev/catmulti-onesession.R.
# A package's whole suite in ONE R session, alphabetical, as covr and
# R CMD check run it. Args: package, arm (lane or base; base runs the
# worktree tests of the package on rellib-r9).
# Output: dev/catmulti-rev-out/onesession-<pkg>[-base].txt
p <- commandArgs(TRUE)[1]
arm <- commandArgs(TRUE)[2]
if (is.na(arm)) arm <- "lane"
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-catmulti-lib", "C:/Users/adf44/source/r/rellib-r9", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win", NOT_CRAN = "true",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-catmulti/dev/stan-cache")
suppressPackageStartupMessages({library(testthat); library(p, character.only = TRUE)})
cat("lib:", dirname(find.package("frmtmb")), "\n")
d <- if (p == "frmtmb") "tests/testthat" else file.path("extensions", p, "tests/testthat")
r <- as.data.frame(test_dir(d, package = p, load_package = "none", reporter = "silent",
                            stop_on_failure = FALSE))
cat(sprintf("ONESESSION %s pass=%d fail=%d err=%d skip=%d warn=%d files=%d\n", p,
            sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped), sum(r$warning),
            length(unique(r$file))))
bad <- r[r$failed > 0 | r$error, c("file", "test")]
if (nrow(bad)) print(bad)
