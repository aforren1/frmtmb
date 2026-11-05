# Review of lane gpgap, item 8: rerun what the lane's change reaches, on its
# library, one session per suite, in file order: the gp files of core and the
# whole spatial suite (whose helper-brms.R the lane changed).
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-rerun.R <log> core|spatial
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(testthat))
if (arm == "core") {
  suppressPackageStartupMessages(library(frmtmb))
  fs <- c("test-brms-lik-gp-exact.R", "test-brms-likelihood.R", "test-brms-agreement.R",
          "test-gp-by.R", "test-gp-multidim.R", "test-gp-nugget-floor.R", "test-v18.R")
  ds <- "tests/testthat"
  pkg <- "frmtmb"
} else {
  suppressPackageStartupMessages(library(frmtmb.spatial))
  ds <- "extensions/frmtmb.spatial/tests/testthat"
  fs <- basename(Sys.glob(file.path(ds, "test-*.R")))
  pkg <- "frmtmb.spatial"
}
cat("core:", dirname(find.package("frmtmb")), "pkg:", dirname(find.package(pkg)), "\n")
tot <- c(0, 0, 0, 0, 0)
env <- test_env(pkg)
for (f in fs) {
  res <- as.data.frame(test_file(file.path(ds, f), package = pkg, env = env,
                                 reporter = "silent", stop_on_failure = FALSE))
  v <- c(sum(res$passed), sum(res$failed), sum(as.integer(res$error)),
         sum(as.integer(res$skipped)), sum(res$warning))
  tot <- tot + v
  cat(sprintf("RERUN %s pass=%d fail=%d error=%d skip=%d warn=%d\n", f, v[1], v[2], v[3], v[4],
              v[5]))
}
cat(sprintf("RERUN total pass=%d fail=%d error=%d skip=%d warn=%d\n", tot[1], tot[2], tot[3],
            tot[4], tot[5]))
