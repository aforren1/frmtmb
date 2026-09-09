# lane tmbstan, probe 03: the exact inline skip lines to rewrite, with
# the test_that block each belongs to and whether that block samples.
dir <- paste0("C:/Users/adf44/source/r/frmtmb-wt-tmbstan/",
              "extensions/frmtmb.sample/tests/testthat")
files <- c("test-sample-direct.R", "test-sampling-ported.R",
           "test-parallel-chains.R", "test-scale.R",
           "test-draws-methods.R", "test-loo.R", "test-reparam.R")
for (f in files) {
  ln <- readLines(file.path(dir, f))
  starts <- grep("^test_that[(]", ln)
  cat("\n#### ", f, "  (", length(ln), " lines)\n", sep = "")
  hits <- grep('skip_if_not_installed[(]"(tmbstan|rstan)"[)]', ln)
  for (h in hits) {
    b <- suppressWarnings(max(starts[starts < h]))
    cat(sprintf("  L%-4d block L%-4d | %s\n", h,
                if (is.finite(b)) b else 0L, trimws(ln[h])))
  }
}
