## The PAIRED statistic, measured on the shipped generator.
##
## dev/coh-absent.R killed the first candidate. Comparing the mean of
## four treatment ratios against the largest of four null ratios PASSED
## on data whose truth had no subject-by-condition effect at all
## (1.0889 against 1.0772), because two independent blocks of four
## drawn from the same distribution land either way about a quarter of
## the time. An assertion that passes when the thing it guards is
## absent is asserting nothing.
##
## The pairing is available and costs nothing. The generator draws the
## id-by-condition deviations at sd 1 and scales, so at one seed the
## treatment and null data differ ONLY by that scaling: with the effect
## absent from both, the two arms are the same data and the per-seed
## quotient is exactly 1. So `min(treat / null) > 1` fails
## deterministically where the mean-against-max statistic flipped a
## coin.
##
## This measures the quotient over 40 seeds, which is what says whether
## the four the test pins are typical.
##
## Run: COH_SEEDS=40 Rscript dev/coh-paired.R
.libPaths(c("C:/Users/adf44/source/r/coh-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.coupling))
setwd("extensions/frmtmb.coupling/tests/testthat")
test_that <- function(...) invisible(NULL)
source("test-coherence.R", local = TRUE)

n <- as.integer(Sys.getenv("COH_SEEDS", "40"))
out <- file.path("..", "..", "..", "..", "dev", "coh-paired.tsv")
seeds <- 2600L + seq_len(n)
one <- function(sd_ic, s) {
  d <- cp_cells(s, n_sub = 16L, n_rep = 6L, sd_idcond = sd_ic)
  cp_width("cond + (1 | id) + (1 | id:cond)", d) /
    cp_width("cond + (1 | id)", d)
}
tr <- numeric(n); nu <- numeric(n)
for (i in seq_len(n)) {
  tr[i] <- one(0.5, seeds[i])
  nu[i] <- one(0, seeds[i])
  cat(paste(c(paste0("seed=", seeds[i]),
              paste0("treat=", signif(tr[i], 8)),
              paste0("null=", signif(nu[i], 8)),
              paste0("quotient=", signif(tr[i] / nu[i], 8))),
            collapse = "\t"), "\n", sep = "", file = out, append = TRUE)
}
q <- tr / nu
cat(sprintf("\n%d seeds. treat [%.3f, %.3f]  null [%.3f, %.3f]\n", n,
            min(tr), max(tr), min(nu), max(nu)))
cat(sprintf(paste0("quotient min %.3f  median %.3f  max %.3f;",
                   "  seeds at or below 1: %d\n"),
            min(q), stats::median(q), max(q), sum(q <= 1)))
cat(sprintf("the four seeds the test pins, 2610:2613: %s\n",
            paste(signif(q[seeds %in% 2610:2613], 4), collapse = " ")))
