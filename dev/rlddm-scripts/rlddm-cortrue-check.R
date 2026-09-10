# PUNCH ROUND 2. The test file's own cor_true, at the tier's size,
# against dev/rlddm-scripts/rlddm-cortrue.R's independent computation.
#
#   Rscript dev/rlddm-scripts/rlddm-cortrue-check.R <lib>
#
# Seed 20260908, the tier's own design. NO FIT: this runs the test
# file's `learn_rlddm_data()` and `learn_cor_max()` at 100 by 200 and
# compares the number the row will record against the number computed
# from scratch by a script that shares no code with it.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(testthat)
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.learn)
})
# The tier is GATED OFF while the file is sourced: sourcing it runs its
# test_that() blocks, and with the gate open that is a 20,000-row fit
# this check does not need. FRMTMB_SCALE_ROW naming no row makes both
# rows skip while the function definitions still land.
Sys.setenv(FRMTMB_SCALE_TESTS = "true",
           FRMTMB_SCALE_ROW = "no-such-row")

here <- "extensions/frmtmb.learn/tests/testthat"
e <- new.env(parent = globalenv())
sys.source(file.path(here, "helper-scale.R"), envir = e)
sys.source(file.path(here, "test-scale.R"), envir = e)

d <- e$learn_rlddm_data()
cat("rows:", nrow(d), " learners:", nlevels(d$id), "\n")
cat("cor_true the row will record:",
    sprintf("%.6f", e$learn_cor_max(attr(d, "dev_fitted_par"))), "\n")
cat("the same truths as DRAWN    :",
    sprintf("%.6f", e$learn_cor_max(attr(d, "dev_drawn"))), "\n")
cat("mean(own floor)             :",
    format(mean(attr(d, "own_floor")), digits = 9), "\n")
