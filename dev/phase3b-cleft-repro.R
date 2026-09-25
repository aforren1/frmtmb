# Punch round 1: the cleft recovery arm stopped at "NA/NaN gradient
# evaluation". Fit the cleft design at 8 subjects by 400 trials over
# seeds, on the library named first, and report which fail.
# Usage: Rscript dev/phase3b-cleft-repro.R <lib> <seed_from> <seed_to>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
for (seed in seq(as.integer(a[2]), as.integer(a[3]))) {
  set.seed(seed)
  NS <- 8L; NT <- 400L
  u <- rnorm(NS, 0, 0.35); b <- rnorm(NS, 0, 0.2)
  d <- do.call(rbind, lapply(seq_len(NS), function(s) {
    cond <- rep(0:1, length.out = NT)
    x <- ddm_simulate(NT, mu = 0.4 + 0.9 * cond + u[s],
                      bs = 1.4 * exp(b[s]), ndt = 0.25, bias = 0.5)
    x$cond <- cond; x$s <- factor(s); x
  }))
  d$code <- 0L; d$y2 <- d$rt
  fast <- d$rt < 0.45
  d$code[fast] <- -1L; d$rt[fast] <- 0.45
  pick <- which(!fast)[seq(5, sum(!fast), by = 5)]
  lo <- floor(d$rt[pick] * 10) / 10
  d$code[pick] <- 2L; d$y2[pick] <- lo + 0.1; d$rt[pick] <- pmax(lo, 0.45)
  r <- tryCatch({
    fit <- frm(bf(rt | dec(upper) + cens(code, y2) ~ cond + (1 | s),
                  bs ~ 1 + (1 | s), ndt ~ 1, bias = 0.5),
               family = wiener(), data = d)
    sprintf("ok code %d logLik %.10f", fit$opt$convergence, as.numeric(logLik(fit)))
  }, error = function(e) substr(conditionMessage(e), 1, 60))
  cat(seed, r, "\n")
}
