# Lane arstep: does the scan need its AR grid (0, 0.5, 0.9, 1)? The `ar()`
# step rows of the seed scan and 30 more `ar()` rows (curves step and mix,
# n 40, 80, 160, seeds 1 to 5; data as dev/review-ncgp-r3-ng.R), fitted with
# the grid and with 0 only (a namespace mutant), in one process, arms
# interleaved row by row. Args: arm
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
rev_libs(args[1])
ns <- asNamespace("frmtmb")
orig <- ns$gp_marginal_opt
txt <- paste(deparse(orig), collapse = "\n")
t2 <- sub("c(0, 0.5, 0.9, 1)", "c(0)", txt, fixed = TRUE)
stopifnot(t2 != txt)
zero <- eval(parse(text = t2)); environment(zero) <- ns
put <- function(f) { unlockBinding("gp_marginal_opt", ns); assign("gp_marginal_opt", f, envir = ns); lockBinding("gp_marginal_opt", ns) }
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x, step = function(x) 2 * (x > 5) - 1)
rows <- list(c("step", 160, 2), c("step", 80, 22), c("step", 80, 23), c("step", 80, 4), c("step", 40, 2))
for (cn in names(curves)) for (n in c(40, 80, 160)) for (seed in 1:5) rows[[length(rows) + 1]] <- c(cn, n, seed)
for (rw in rows) {
  cn <- rw[1]; n <- as.integer(rw[2]); seed <- as.integer(rw[3])
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  d <- data.frame(x = x, t = seq_len(n), y = curves[[cn]](x) + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03)))
  out <- sprintf("PG %s %d %d", cn, n, seed)
  for (arm in c("grid", "zero")) {
    put(if (arm == "grid") orig else zero)
    t0 <- proc.time()[[3]]
    r <- rev_fit(y ~ gp(x) + ar(time = t), d, family = gaussian())
    el <- proc.time()[[3]] - t0
    out <- paste(out, sprintf("| %s ll %.6f conv %d w %d wall %.2f", arm,
                              if (inherits(r$fit, "error")) NA else as.numeric(logLik(r$fit)),
                              if (inherits(r$fit, "error")) NA else r$fit$opt$convergence, length(r$warn), el))
  }
  cat(out, "\n"); flush.console()
}
put(orig)
cat("DONE\n")
