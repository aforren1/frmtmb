# Backlog 4002: does covr install frmtmb.spatial from inside an R process
# on this machine? Five installs alone (type = "none"), then one run of the
# package's tests under covr. R_LIBS (set by the caller) holds the lane's
# library, rellib-r21 and the user library, so covr's child install finds
# frmtmb.
.libPaths(strsplit(Sys.getenv("R_LIBS"), ";")[[1]])
Sys.setenv(NOT_CRAN = "true")
cat("covr", format(packageVersion("covr")), "R", R.version.string, "\n")
cat(".libPaths():", .libPaths(), sep = "\n  ")
pkg <- "C:/Users/adf44/source/r/frmtmb-wt-arspat2/extensions/frmtmb.spatial"
ok <- 0
for (i in 1:5) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch({
    covr::package_coverage(pkg, type = "none", quiet = TRUE)
    "installed"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok <- ok + identical(r, "installed")
  cat(sprintf("attempt %d (type none): %s, %.1f s\n", i, substr(r, 1, 300),
              proc.time()[["elapsed"]] - t0))
}
cat(sprintf("installs that worked: %d of 5\n", ok))
t0 <- proc.time()[["elapsed"]]
cv <- tryCatch(covr::package_coverage(pkg, type = "tests", quiet = TRUE),
               error = function(e) e)
if (inherits(cv, "error")) {
  cat("type tests: ERROR:", conditionMessage(cv), "\n")
} else {
  cat(sprintf("type tests: coverage %.2f percent, %.1f s\n",
              covr::percent_coverage(cv), proc.time()[["elapsed"]] - t0))
}
