# Reviewer: CPU load for dev/nanse-rev2-determinism.R. Busy-loops for
# the given seconds.
#   Rscript dev/nanse-rev2-burn.R <seconds>
s <- as.numeric(commandArgs(TRUE)[1])
t0 <- proc.time()[["elapsed"]]
x <- 0
while (proc.time()[["elapsed"]] - t0 < s) x <- x + sum(sqrt(seq_len(1e5)))
