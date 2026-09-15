## Re-check spot: does dev/coh-absent.R's mutation probe mutate what it
## claims to? A gsub over a deparsed body that FAILS to match leaves
## the generator untouched, and the probe would then be comparing the
## shipped code with itself. This checks the mutated body really
## carries the other spelling before believing the FALSE it prints.
.libPaths(c("C:/Users/adf44/source/r/cohrev-lib",
            "C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
setwd("extensions/frmtmb.coupling/tests/testthat")
test_that <- function(...) invisible(NULL)
source("test-coherence.R", local = TRUE)

mutated <- cp_cells
body(mutated) <- parse(text = gsub(
  "stats::rnorm(n_sub * 2L) * sd_idcond",
  "stats::rnorm(n_sub * 2L, 0, sd_idcond)",
  paste(deparse(body(cp_cells)), collapse = "\n"), fixed = TRUE))[[1L]]

src_ship <- paste(deparse(body(cp_cells)), collapse = "\n")
src_mut <- paste(deparse(body(mutated)), collapse = "\n")
cat("the bodies differ at all:              ",
    !identical(src_ship, src_mut), "[wanted TRUE]\n")
cat("shipped body has `rnorm(n_sub * 2L) * sd_idcond`:",
    grepl("stats::rnorm(n_sub * 2L) * sd_idcond", src_ship,
          fixed = TRUE), "[wanted TRUE]\n")
cat("mutated body has `rnorm(n_sub * 2L, 0, sd_idcond)`:",
    grepl("stats::rnorm(n_sub * 2L, 0, sd_idcond)", src_mut,
          fixed = TRUE), "[wanted TRUE]\n")
cat("mutated body still has the shipped spelling:",
    grepl("stats::rnorm(n_sub * 2L) * sd_idcond", src_mut,
          fixed = TRUE), "[wanted FALSE]\n")
cat("only that one line changed:",
    sum(strsplit(src_ship, "\n")[[1L]] !=
          strsplit(src_mut, "\n")[[1L]]), "[wanted 1]\n")

stream_after <- function(gen, sd_ic) {
  gen(2610L, n_sub = 4L, n_rep = 2L, sd_idcond = sd_ic)
  .Random.seed
}
cat("\nshipped, same stream at sd 0 and 0.5:",
    identical(stream_after(cp_cells, 0), stream_after(cp_cells, 0.5)),
    "[wanted TRUE]\n")
cat("mutated, same stream at sd 0 and 0.5:",
    identical(stream_after(mutated, 0), stream_after(mutated, 0.5)),
    "[wanted FALSE]\n")

## And the consequence the probe exists to demonstrate: with the
## mutation in place the two arms are no longer the same data, so the
## absent case would no longer give a quotient of exactly 1.
a <- mutated(2610L, n_sub = 16L, n_rep = 6L, sd_idcond = 0)
b <- cp_cells(2610L, n_sub = 16L, n_rep = 6L, sd_idcond = 0)
cat("\nmutated and shipped generators agree at sd 0:",
    identical(a, b), "[wanted FALSE: the mutation moves the draws]\n")
