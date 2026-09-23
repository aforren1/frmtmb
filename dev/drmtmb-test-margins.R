# Measures, for every expect_same_model() call in
# tests/testthat/test-drmtmb-agreement.R, the quantities it bounds, so the
# margins quoted in that file's header come from its own models and
# seeds. The test file is parsed; its helpers are evaluated as written and
# each test body runs with expect_same_model() swapped for a recorder.
.libPaths(c("C:/Users/adf44/source/r/drmtmb-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(testthat); library(frmtmb) })
ex <- parse(file.path("C:/Users/adf44/source/r/frmtmb-wt-drmtmb",
                      "tests/testthat/test-drmtmb-agreement.R"))
env <- new.env(parent = globalenv())
is_test <- function(e) is.call(e) && identical(e[[1]], as.name("test_that"))
for (e in ex) if (!is_test(e)) eval(e, env)
rec <- list()
for (e in ex) if (is_test(e)) {
  label <- e[[2]]
  tenv <- new.env(parent = env)
  tenv$skip_unless_drmtmb <- function() invisible()
  k <- 0
  tenv$expect_same_model <- function(fd, ff, map) {
    k <<- k + 1
    m <- env$drm_measure(fd, ff, map)
    s <- abs(m$ll_frm)
    rec[[length(rec) + 1]] <<- data.frame(
      test = substr(label, 1, 40), call = k,
      ll = abs(m$ll_frm - m$ll_drm) / s, at_drm = abs(m$gap_at_drm) / s,
      off = abs(m$gap_off) / s, diff_se = m$diff_over_se,
      log_se = m$log_se_ratio)
  }
  eval(e[[3]], tenv)
}
r <- do.call(rbind, rec)
op <- options(width = 160); print(r, digits = 3); options(op)
# One line per expectation in expect_same_model(), so the header of the
# test file can quote all five margins and not a subset.
line <- function(label, worst, bound) {
  cat(sprintf("%-34s worst %8.2e  bound %8.1e  margin %5.1f x\n",
              label, worst, bound, bound / worst))
}
cat("\n")
line("gap between the two optima", max(r$ll), 1e-9)
line("gap at drmTMB's optimum", max(r$at_drm), 1e-9)
line("gap one SE off the optimum", max(r$off), 1e-11)
line("estimate gap / standard error", max(r$diff_se), 1e-2)
line("|log standard-error ratio|", max(r$log_se), 1e-3)
