## The case where the guarded thing is ABSENT.
##
## The shipped assertion is not a bug pin, so there is no unfixed code
## to run it against. What it can be run against is data whose truth
## has no subject-by-condition effect at all: the statistic must then
## FAIL, or it is asserting nothing. The mock treatment arm below is
## drawn at sd(id:cond) = 0 on four seeds the null arm does not use, so
## the two arms are independent draws of the same null truth.
##
## Run: Rscript dev/coh-absent.R
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
setwd("extensions/frmtmb.coupling/tests/testthat")
test_that <- function(...) invisible(NULL)
source("test-coherence.R", local = TRUE)

ratio <- function(sd_ic, seeds) {
  vapply(seeds, function(s) {
    d <- cp_cells(s, n_sub = 16L, n_rep = 6L, sd_idcond = sd_ic)
    cp_width("cond + (1 | id) + (1 | id:cond)", d) /
      cp_width("cond + (1 | id)", d)
  }, numeric(1))
}
seeds <- 2610:2613
mock <- ratio(0, 2620:2623)
null <- ratio(0, seeds)
cat("mock treatment (sd_idcond = 0):",
    paste(signif(mock, 5), collapse = " "), "\n")
cat("null           (sd_idcond = 0):",
    paste(signif(null, 5), collapse = " "), "\n")

## 1. The statistic that was NOT shipped, on independent blocks of the
##    same null truth. It passes, which is why it was not shipped.
cat("\nunpaired, mean(mock) =", mean(mock), " max(null) =", max(null),
    "\n")
cat("  mean(treat) > max(null) [wanted FALSE]:",
    mean(mock) > max(null), "\n")

## 2. The statistic that ships, at the SAME seeds, which is what pairing
##    means. With the effect absent from both arms the two calls draw
##    the same data, so the quotient is exactly 1 and the assertion
##    fails by construction rather than by luck.
same <- ratio(0, seeds)
cat("\npaired at the same seeds, quotient:",
    paste(signif(same / null, 17), collapse = " "), "\n")
cat("  identical widths:", identical(same, null), "\n")
cat("  min(treat / null) > 1 [wanted FALSE]:", min(same / null) > 1,
    "\n")
cat("  min(treat) > 1 [wanted FALSE]:", min(same) > 1, "\n")

## 3. The guard's own guard: the stream assertion the test now carries,
##    run against the ONE spelling that would break the pairing. No fit.
stream_after <- function(gen, sd_ic) {
  gen(2610L, n_sub = 4L, n_rep = 2L, sd_idcond = sd_ic)
  .Random.seed
}
mutated <- cp_cells
body(mutated) <- parse(text = gsub(
  "stats::rnorm(n_sub * 2L) * sd_idcond",
  "stats::rnorm(n_sub * 2L, 0, sd_idcond)",
  paste(deparse(body(cp_cells)), collapse = "\n"),
  fixed = TRUE))[[1L]]
cat("\nshipped cp_cells, same stream at sd 0 and 0.5:",
    identical(stream_after(cp_cells, 0), stream_after(cp_cells, 0.5)),
    "[wanted TRUE]\n")
cat("mutated to rnorm(n, 0, sd), same stream:",
    identical(stream_after(mutated, 0), stream_after(mutated, 0.5)),
    "[wanted FALSE]\n")
