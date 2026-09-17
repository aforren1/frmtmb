# Reviewer, lane wt-conditions (round 2): the lane's census test on each
# of the 8 BASE source trees (main checkout, at e031e8c). Every one must
# fail.
#   Rscript dev/conditions-rev-census-base.R
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-conditions"
main <- "C:/Users/adf44/source/r/frmtmb"
scratch <- file.path(wt, "dev/conditions-rev-census-base")
unlink(scratch, recursive = TRUE); dir.create(scratch)
cat("R/windows or R/unix dirs in the lane trees:",
    length(Sys.glob(file.path(wt, c("R", "extensions/*/R"),
                              c("windows", "unix")))), "\n")
cat("R/ files not ending in .R in the lane trees:",
    length(grep("[.]R$", list.files(Sys.glob(file.path(
      wt, c("R", "extensions/*/R")))), invert = TRUE, value = TRUE)), "\n")
pk <- c(frmtmb = "", setNames(paste0("extensions/", list.files(
  file.path(wt, "extensions"))), list.files(file.path(wt, "extensions"))))
for (p in names(pk)) {
  d <- file.path(scratch, p)
  dir.create(file.path(d, "tests", "testthat"), recursive = TRUE)
  src <- file.path(main, pk[[p]])
  file.copy(file.path(src, "R"), d, recursive = TRUE)
  if (dir.exists(file.path(src, "inst")))
    file.copy(file.path(src, "inst"), d, recursive = TRUE)
  file.copy(file.path(src, "DESCRIPTION"), d)
  file.copy(file.path(wt, pk[[p]], "tests/testthat/test-conditions-census.R"),
            file.path(d, "tests", "testthat"))
  r <- as.data.frame(test_file(file.path(d, "tests/testthat/test-conditions-census.R"),
                               reporter = "silent", package = "frmtmb"))
  cat(sprintf("%-16s pass %d fail %d error %d skip %d\n", p, sum(r$passed),
              sum(r$failed), sum(r$error), sum(r$skipped)))
}
unlink(scratch, recursive = TRUE)
