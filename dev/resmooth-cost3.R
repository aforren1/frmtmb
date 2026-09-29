# Lane wt-resmooth, nits round. The reviewer's four cost cells plus the
# one the batching fix is for, on three arms: 0.64.0, this lane before
# the nits round, and after it.
#
# The clock instrument is the reviewer's (proc.time() ticks at 10 ms
# here, so each arm grows a block of repeats past 1.2 s and reports the
# minimum of three blocks over the repeat count) and the CONTROL is
# carried, because the arms are separate processes and cannot be
# interleaved across two installed versions of one package. The clock is
# the WEAKER instrument here and is reported beside the one that does not
# move with load: the exact count of model evaluations the
# finite-difference route makes, which dev/resmooth-batchcost.R reports
# for the two lane arms.
#   COST3=base   Rscript dev/resmooth-cost3.R > dev/resmooth-cost3-base.txt
#   COST3=before Rscript dev/resmooth-cost3.R > dev/resmooth-cost3-before.txt
#   COST3=after  Rscript dev/resmooth-cost3.R > dev/resmooth-cost3-after.txt
arm <- Sys.getenv("COST3")
lib <- switch(arm,
              base = "C:/Users/adf44/source/r/rellib-r3",
              before = "C:/Users/adf44/source/r/wt-resmooth-lib",
              after = "C:/Users/adf44/source/r/wt-resmooth-lib2",
              stop("set COST3 to base, before or after"))
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM:", arm, "from", find.package("frmtmb"), "\n")

block <- function(f, target = 1.2, rounds = 3) {
  reps <- 1L
  repeat {
    t0 <- proc.time()
    for (i in seq_len(reps)) f()
    e <- (proc.time() - t0)[["elapsed"]]
    if (e >= target) break
    reps <- reps * 2L
    if (reps > 4096L) break
  }
  best <- Inf
  for (r in seq_len(rounds)) {
    t0 <- proc.time()
    for (i in seq_len(reps)) f()
    best <- min(best, (proc.time() - t0)[["elapsed"]])
  }
  best / reps
}
say <- function(lab, f) cat(sprintf("%-46s %.4f s\n", lab, block(f)))

set.seed(41)
nD <- 160
dD <- data.frame(x = stats::runif(nD, -2, 2))
latD <- 1.1 * sin(2 * dD$x) + stats::rlogis(nD)
dD$y <- factor(cut(latD, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fD <- suppressWarnings(frm(bf(y ~ gp(x)), family = cumulative(), data = dD))
ndD <- data.frame(x = c(-1.5, 0, 1.5))
say("A gp(x) n_b=160, newdata 3 NEW positions, NA",
    function() fitted(fD, newdata = ndD, re_formula = NA))
say("B gp(x) n_b=160, in sample, NA",
    function() fitted(fD, re_formula = NA))

set.seed(31)
ng <- 40
per <- 20
dB <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                 x = stats::runif(ng * per, -2, 2))
latB <- sin(1.5 * dB$x) + stats::rnorm(ng, 0, 0.8)[dB$g] +
  stats::rlogis(nrow(dB))
dB$y <- factor(cut(latB, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fB <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | g)),
                           family = cumulative(), data = dB))
say("D s(x,k=8)+(1|g) 40 lv, in sample, NULL",
    function() fitted(fB, re_formula = NULL))
say("E s(x,k=8)+(1|g) 40 lv, in sample, NA",
    function() fitted(fB, re_formula = NA))

set.seed(59)
dH <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                 x = stats::runif(ng * per, -2, 2))
latH <- 0.8 * dH$x + stats::rnorm(ng, 0, 0.7)[dH$g] +
  stats::rlogis(nrow(dH))
dH$y <- factor(cut(latH, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fH <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = cumulative(),
                           data = dH))
say("CONTROL x+(1|g) 40 lv, in sample, NULL",
    function() fitted(fH, re_formula = NULL))
cat("DONE\n")
