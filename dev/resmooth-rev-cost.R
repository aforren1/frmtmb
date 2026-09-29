# Reviewer: does differencing every smooth's coefficients cost wall clock?
# proc.time() ticks at 10 ms here, so each arm runs a BLOCK of repeats
# until it passes 1.2 s, and the reported figure is the minimum of three
# blocks divided by the repeat count. A CONTROL (a no-smooth ordinal fit,
# which this change cannot touch) must report a ratio near 1.
ARM <- if (identical(Sys.getenv("REVLIB"), "base")) "base" else "lane"
LIB <- if (ARM == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM:", ARM, "\n")

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
  c(per_call = best / reps, reps = reps)
}

set.seed(41)
nD <- 160
dD <- data.frame(x = stats::runif(nD, -2, 2))
latD <- 1.1 * sin(2 * dD$x) + stats::rlogis(nD)
dD$y <- factor(cut(latD, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fD <- suppressWarnings(frm(bf(y ~ gp(x)), family = cumulative(), data = dD))
ndD <- data.frame(x = c(-1.5, 0, 1.5))
r <- block(function() fitted(fD, newdata = ndD, re_formula = NA))
cat(sprintf("gp(x), n_b = %d : %.4f s per fitted() (reps %d)\n",
            length(fD$estimates[["b"]]), r[["per_call"]], r[["reps"]]))

set.seed(31)
ng <- 40; per <- 20
dB <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                 x = stats::runif(ng * per, -2, 2))
latB <- sin(1.5 * dB$x) + stats::rnorm(ng, 0, 0.8)[dB$g] +
  stats::rlogis(nrow(dB))
dB$y <- factor(cut(latB, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fB <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | g)),
                           family = cumulative(), data = dB))
r <- block(function() fitted(fB, re_formula = NULL))
cat(sprintf("s(x,k=8)+(1|g) 40 levels in sample NULL : %.4f s (reps %d)\n",
            r[["per_call"]], r[["reps"]]))
r <- block(function() fitted(fB, re_formula = NA))
cat(sprintf("s(x,k=8)+(1|g) 40 levels in sample NA   : %.4f s (reps %d)\n",
            r[["per_call"]], r[["reps"]]))

# CONTROL: no smooth at all, so neither arm differences anything new
set.seed(59)
dH <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                 x = stats::runif(ng * per, -2, 2))
latH <- 0.8 * dH$x + stats::rnorm(ng, 0, 0.7)[dH$g] +
  stats::rlogis(nrow(dH))
dH$y <- factor(cut(latH, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fH <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = cumulative(),
                           data = dH))
r <- block(function() fitted(fH, re_formula = NULL))
cat(sprintf("CONTROL x+(1|g) 40 levels in sample NULL: %.4f s (reps %d)\n",
            r[["per_call"]], r[["reps"]]))
cat("DONE\n")
