# The tests of fix round 2 read against the base build. A test that pins a
# defect of the base must fail there; one that pins a defect this lane brought
# in cannot, and the report says which is which. Log impsamp-out/r16check.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
})
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
tests <- c("a threshold does not step by its own size",
           "a sd far from one reads its own precision",
           "the rounds start where the Laplace route ends",
           "a predictor in small response units steps in those units",
           "a count whose worst group is below the floor is not named")
for (d in tests) {
  r <- as.data.frame(test_file("tests/testthat/test-importance.R", desc = d,
                               reporter = "silent", stop_on_failure = FALSE,
                               package = "frmtmb"))
  cat(sprintf("%-60s pass %d fail %d error %d\n", d, sum(r$passed),
              sum(r$failed), sum(r$error)))
}
