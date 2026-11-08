# Lane arstep: what the AR marginal costs at size. y ~ gp(x) + ar(time = t)
# on distinct positions, a smooth curve plus AR(0.4) noise sd 0.1; and the
# same without ar() as the control, whose route this lane does not change.
# frm() wall time and ll, `reps` rounds, interleaved by size.
# Seed 9300 + n. Args: arm sizes(comma) reps
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
arm <- args[1]; sizes <- as.integer(strsplit(args[2], ",")[[1]]); reps <- as.integer(args[3])
rev_libs(arm)
for (rep in seq_len(reps)) for (n in sizes) {
  set.seed(9300 + n)
  x <- sort(runif(n, 0, 10))
  d <- data.frame(x = x, t = seq_len(n),
                  y = sin(x) + 0.5 * sin(3 * x) + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.1)))
  for (what in c("ar", "control")) {
    f <- if (what == "ar") y ~ gp(x) + ar(time = t) else y ~ gp(x)
    gc()
    t0 <- proc.time()[[3]]
    r <- rev_fit(f, d, family = gaussian())
    el <- proc.time()[[3]] - t0
    if (inherits(r$fit, "error")) { cat(sprintf("ACO %s %s n %d rep %d ERROR %s\n", arm, what, n, rep, conditionMessage(r$fit))); next }
    cat(sprintf("ACO %s %s n %d rep %d wall %.2f ll %.6f conv %d w %d marginal %s\n", arm, what, n, rep, el,
                as.numeric(logLik(r$fit)), r$fit$opt$convergence, length(r$warn),
                !is.null(r$fit$obj$env$.frmtmb_gp_marginal)))
    flush.console()
  }
}
cat("DONE\n")
