## Item 4 of the review: construct the ABSENT case against the SHIPPED
## assertion, independently of dev/coh-absent.R.
##
## The lane's first candidate passed with the effect absent. Its
## replacement is `expect_gt(min(treat / null), 1)`. This runs that
## expression, verbatim, with the truth carrying NO subject-by-
## condition effect in either arm, and it runs it through testthat so
## that a failure is a failure and not a printed FALSE.
##
## Three constructions, because one seed set is not a construction:
##   A  the test's own seeds, 2610:2613, effect absent
##   B  four seeds the lane never used, 2700:2703, effect absent
##   C  the effect present but small, sd(id:cond) = 0.15, to find where
##      the guard stops discriminating
.libPaths(c("C:/Users/adf44/source/r/cohrev-lib",
            "C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
suppressMessages(library(testthat))
setwd("extensions/frmtmb.coupling/tests/testthat")
test_that <- function(...) invisible(NULL)
source("test-coherence.R", local = TRUE)
rm(test_that)

ratio <- function(sd_ic, seeds) {
  vapply(seeds, function(s) {
    d <- cp_cells(s, n_sub = 16L, n_rep = 6L, sd_idcond = sd_ic)
    cp_width("cond + (1 | id) + (1 | id:cond)", d) /
      cp_width("cond + (1 | id)", d)
  }, numeric(1))
}

## The shipped assertion, lifted out of the test file so that the same
## expression is used in every arm below.
shipped <- function(treat, null, label) {
  r <- testthat::test_that(label, {
    testthat::expect_gt(min(treat / null), 1)
    testthat::expect_gt(min(treat), 1)
  })
  cat(sprintf("%-40s PASSED: %s\n", label, isTRUE(r)))
  invisible(r)
}

report <- function(tag, treat, null) {
  cat(sprintf("\n%s\n  treat: %s\n  null : %s\n  quotient: %s\n", tag,
              paste(signif(treat, 6), collapse = " "),
              paste(signif(null, 6), collapse = " "),
              paste(signif(treat / null, 8), collapse = " ")))
  cat(sprintf("  widths bitwise identical between arms: %s\n",
              identical(treat, null)))
  cat(sprintf("  min(treat/null) = %.17g;  min(treat) = %.17g\n",
              min(treat / null), min(treat)))
}

cat("==== A: the test's own seeds, effect ABSENT in both arms ====\n")
nullA <- ratio(0, 2610:2613)
treatA <- ratio(0, 2610:2613)
report("A", treatA, nullA)
shipped(treatA, nullA, "A absent, test's own seeds [wanted FAIL]")

cat("\n==== B: four unused seeds, effect ABSENT in both arms ====\n")
nullB <- ratio(0, 2700:2703)
treatB <- ratio(0, 2700:2703)
report("B", treatB, nullB)
shipped(treatB, nullB, "B absent, fresh seeds [wanted FAIL]")

cat("\n==== C: effect PRESENT but small, sd(id:cond) = 0.15 ====\n")
treatC <- ratio(0.15, 2700:2703)
report("C", treatC, nullB)
shipped(treatC, nullB, "C small effect [informational]")

cat("\n==== D: the shipped arm, sd(id:cond) = 0.5, own seeds ====\n")
treatD <- ratio(0.5, 2610:2613)
report("D", treatD, nullA)
shipped(treatD, nullA, "D shipped arm [wanted PASS]")
