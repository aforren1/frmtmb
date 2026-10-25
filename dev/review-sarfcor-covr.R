# Reviewer: one covr::package_coverage() attempt on frmtmb.spatial,
# install_path set, as test-coverage.yaml calls it.
# Rscript dev/review-sarfcor-covr.R base|lane
arm <- commandArgs(TRUE)[1]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor"
libs <- c(if (arm == "lane") "C:/Users/adf44/source/r/wt-sarfcor-lib",
          "C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true")
src <- if (arm == "lane") file.path(wt, "extensions/frmtmb.spatial") else
  file.path(wt, "dev/review-sarfcor-covr-base/extensions/frmtmb.spatial")
ip <- file.path(wt, paste0("dev/review-sarfcor-covr-ip-", arm))
unlink(ip, recursive = TRUE)
cat("arm", arm, "covr", format(packageVersion("covr")), "frmtmb",
    find.package("frmtmb"), "\n")
t0 <- proc.time()[["elapsed"]]
cov <- tryCatch(covr::package_coverage(src, quiet = FALSE, clean = FALSE,
                                       install_path = ip),
                error = function(e) e)
cat("elapsed", proc.time()[["elapsed"]] - t0, "\n")
if (inherits(cov, "error")) {
  cat("COVR ERROR:", conditionMessage(cov), "\n")
} else {
  cat("COVR OK: coverage", format(covr::percent_coverage(cov), digits = 4), "%\n")
}
