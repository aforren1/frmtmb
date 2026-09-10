# Lane `latent`, item 2.3: what does `hmm_starts()` cost at the
# realistic scale?
#
# The design is the plan's own: 50 sequences x 500 steps, K = 3
# gaussian, `tr12 ~ (1 | id)`, 25 000 rows. A refit is a whole fit, so
# the answer is n times one fit, and the point of measuring rather than
# asserting that is to say what one fit is HERE, today, on this machine.
#
# HOW THE MEASUREMENT IS INSTRUMENTED. proc.time() ticks at 10.0 ms on
# this machine and two timing claims in an earlier round evaporated
# under replication, so:
#
#  - the blocks are minutes long, far past any clock resolution;
#  - the SAME work is run three times, from the same seed, in ONE
#    process, and the ratio of the slowest block to the fastest is the
#    CONTROL. Identical work through an identical path must report 1.0,
#    and a control far from 1.0 says the machine and not the model is
#    what the row measured;
#  - a LOAD-INDEPENDENT count is carried beside the clock. The same
#    seed gives the same starting values, so the optimizer takes the
#    same path and the gradient-evaluation counts must be IDENTICAL
#    across blocks, not merely close. That is the check the clock
#    cannot make.
#
# Run this in a quiet process: nothing else of this lane's should be
# running.
#
#   Rscript dev/latent-2p3-starts-cost.R [n] [blocks]

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("dev/latent-hmm-sim.R")

args <- commandArgs(trailingOnly = TRUE)
NST <- if (length(args) >= 1L) as.integer(args[[1L]]) else 2L
BLOCKS <- if (length(args) >= 2L) as.integer(args[[2L]]) else 3L
SEED <- 4201L

elapsed <- function(expr) {
  gc(FALSE)
  t0 <- Sys.time()
  force(expr)
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

s <- hmm_sim(seed = 20260910L, ns = 50L, tl = 500L)
d <- s$d
form <- bf(y ~ 1, tr12 ~ 1 + (1 | id))
fam <- hmm(K = 3, gaussian(), time = t, group = id, init = "uniform")
cat("rows:", nrow(d), "  sequences:", nlevels(d$id), "\n")

t_fit <- elapsed(fit <- suppressWarnings(frm(form, family = fam,
                                             data = d)))
cat("one frm() fit, se deferred      :", format(t_fit, digits = 4),
    "s\n")
t_sdr <- elapsed(suppressWarnings(stats::confint(fit)))
cat("the deferred sdreport, once     :", format(t_sdr, digits = 4),
    "s  (hmm_starts() pays this once, for the jitter scale)\n")
cat("outer parameters                :", length(fit$opt$par), "\n")
cat("the fit's own gradient calls    :",
    fit$opt$evaluations[2L], "\n\n")

tt <- numeric(BLOCKS)
gr <- numeric(BLOCKS)
sec_tab <- vector("list", BLOCKS)
for (b in seq_len(BLOCKS)) {
  ms <- NULL
  tt[b] <- elapsed(ms <- hmm_starts(fit, n = NST, jitter = 2,
                                    seed = SEED))
  g <- ms$table$gr_evals[-1L]
  gr[b] <- sum(g, na.rm = TRUE)
  sec_tab[[b]] <- ms$table$seconds[-1L]
  cat(sprintf("block %d: %.2f s for %d refits, %g gradient calls, ",
              b, tt[b], NST, gr[b]))
  cat(sprintf("best %.6f, %d converged, %d not, %d errored\n",
              as.numeric(logLik(ms$best)), ms$n_converged,
              ms$n_not_converged, ms$n_error))
  flush(stdout())
}

cat("\n== the instrument ==\n")
cat("  block seconds        :", format(tt, digits = 5), "\n")
cat("  CONTROL, max / min   :", format(max(tt) / min(tt), digits = 4),
    " (must be near 1.0)\n")
cat("  gradient calls/block :", gr, "\n")
cat("  identical across blocks:",
    length(unique(gr)) == 1L, "\n")

cat("\n== the cost ==\n")
cat("  seconds per refit, minimum over", BLOCKS, "blocks:",
    format(min(tt) / NST, digits = 5), "\n")
cat("  gradient calls per refit       :",
    format(gr[1L] / NST, digits = 6), "\n")
cat("  so hmm_starts(n) costs about   :",
    format(min(tt) / NST, digits = 4), "s x n, plus",
    format(t_sdr, digits = 3), "s once\n")
for (nn in c(4L, 8L, 16L)) {
  cat(sprintf("    n = %2d : %.1f minutes\n", nn,
              (t_sdr + nn * min(tt) / NST) / 60))
}
cat("\nper-refit seconds, every block:\n")
print(round(do.call(rbind, sec_tab), 2))
