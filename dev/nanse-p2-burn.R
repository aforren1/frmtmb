# Lane, punch round 2: CPU load for the determinism and concurrency
# runs. Busy-loops for the given seconds, then exits by itself.
#   Rscript dev/nanse-p2-burn.R <seconds>
s <- as.numeric(commandArgs(TRUE)[1])
t0 <- proc.time()[["elapsed"]]
x <- 0
while (proc.time()[["elapsed"]] - t0 < s) x <- x + sum(sqrt(seq_len(1e5)))
