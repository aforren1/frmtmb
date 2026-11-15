# Review of lane hygiene, item 1252: test-codemeta.R against the 0.82.0
# codemeta.json, the regenerated one, and a tree with no codemeta.json
# (as R CMD check and covr see it). Builds each tree in a temp dir.
#   Rscript dev/reviews/hygiene-review-codemeta.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(testthat))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-hygiene"
arm <- function(label, cm) {
  root <- tempfile("cmtree")
  dir.create(file.path(root, "tests", "testthat"), recursive = TRUE)
  file.copy(file.path(wt, "DESCRIPTION"), root)
  if (!is.null(cm)) file.copy(cm, file.path(root, "codemeta.json"))
  file.copy(file.path(wt, "tests", "testthat", "test-codemeta.R"),
            file.path(root, "tests", "testthat"))
  r <- as.data.frame(test_file(file.path(root, "tests", "testthat",
                                         "test-codemeta.R"),
                               reporter = "silent"))
  cat(sprintf("%-22s pass=%d fail=%d skip=%d\n", label, sum(r$passed),
              sum(r$failed), sum(r$skipped)))
}
arm("0.82.0 codemeta.json", file.path(wt, "dev", "reviews",
                                      "hygiene-review-cm-base.json"))
arm("regenerated", file.path(wt, "codemeta.json"))
arm("no codemeta.json", NULL)
