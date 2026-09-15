## Continuation of dev/cohrev-guard.R: arms B, C and D, with the
## testthat failure caught so that one failing arm does not stop the
## rest, plus two probes the lane did not run.
.libPaths(c("C:/Users/adf44/source/r/cohrev-lib",
            "C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
suppressMessages(library(testthat))

## Probe 1, and it costs nothing: the whole guard rests on the two arms
## drawing the SAME data when the effect is absent, which holds only
## because the generator writes `rnorm(n) * sd` and not
## `rnorm(n, 0, sd)`. Check that R really does skip the draw at sd 0,
## and that the shipped generator therefore consumes the same stream at
## every sd.
set.seed(1); a <- rnorm(3, 0, 0); sa <- .Random.seed
set.seed(1); b <- rnorm(3, 0, 1); sb <- .Random.seed
cat("rnorm(n, 0, 0) SKIPS the draw:", !identical(sa, sb), "\n")
set.seed(1); a2 <- rnorm(3) * 0; sa2 <- .Random.seed
cat("rnorm(n) * 0 SKIPS the draw:  ", !identical(sa2, sb), "\n")

setwd("extensions/frmtmb.coupling/tests/testthat")
test_that <- function(...) invisible(NULL)
source("test-coherence.R", local = TRUE)
rm(test_that)

## Probe 2: does the SHIPPED generator consume the same stream at
## sd_idcond 0 and 0.5? If it ever stops doing so, the two arms are no
## longer the same data and the guard silently becomes the unpaired
## statistic that failed open.
seed_after <- function(sd_ic) {
  cp_cells(2610L, n_sub = 4L, n_rep = 2L, sd_idcond = sd_ic)
  .Random.seed
}
cat("cp_cells consumes the same stream at sd 0 and 0.5:",
    identical(seed_after(0), seed_after(0.5)), "\n")

ratio <- function(sd_ic, seeds) {
  vapply(seeds, function(s) {
    d <- cp_cells(s, n_sub = 16L, n_rep = 6L, sd_idcond = sd_ic)
    cp_width("cond + (1 | id) + (1 | id:cond)", d) /
      cp_width("cond + (1 | id)", d)
  }, numeric(1))
}
shipped <- function(treat, null, label) {
  ok <- tryCatch({
    testthat::test_that(label, {
      testthat::expect_gt(min(treat / null), 1)
    })
    TRUE
  }, error = function(e) FALSE)
  sign_ok <- tryCatch({
    testthat::test_that(label, {
      testthat::expect_gt(min(treat), 1)
    })
    TRUE
  }, error = function(e) FALSE)
  cat(sprintf("%-44s quotient assertion PASSED: %-5s  sign check PASSED: %s\n",
              label, ok, sign_ok))
}
report <- function(tag, treat, null) {
  cat(sprintf("%s\n  treat: %s\n  null : %s\n  quotient: %s\n",
              tag, paste(signif(treat, 6), collapse = " "),
              paste(signif(null, 6), collapse = " "),
              paste(signif(treat / null, 8), collapse = " ")))
  cat(sprintf("  bitwise identical arms: %s;  min(treat/null) = %.17g\n",
              identical(treat, null), min(treat / null)))
}

cat("\n==== B: four unused seeds, effect ABSENT in both arms ====\n")
nullB <- ratio(0, 2700:2703)
treatB <- ratio(0, 2700:2703)
report("B", treatB, nullB)
shipped(treatB, nullB, "B absent, fresh seeds [wanted FAIL]")

cat("\n==== C: effect PRESENT but small, sd(id:cond) = 0.15 ====\n")
treatC <- ratio(0.15, 2700:2703)
report("C", treatC, nullB)
shipped(treatC, nullB, "C small effect [informational]")

cat("\n==== D: the shipped arm, sd(id:cond) = 0.5, fresh seeds ====\n")
treatD <- ratio(0.5, 2700:2703)
report("D", treatD, nullB)
shipped(treatD, nullB, "D shipped arm, fresh seeds [wanted PASS]")
