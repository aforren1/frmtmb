# What the fallback costs in wall clock. Arms are interleaved inside one
# process with a FIXED repeat count per arm (chosen once, in a warmup,
# so no round pays for calibration), the arm order rotates each round,
# each block runs past 1.2 s, and the figure is the MINIMUM over rounds.
#
# Two controls: the same fit entered as a second arm, which must report
# 1.0, and a fixed arithmetic loop that no build changes, whose value
# across builds says how loaded the box was.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")

dd <- make_data(1, FALSE)
m3 <- skew(dd$y)
raw_start <- list(betad = c(log(stats::sd(dd$y)), 2 * sign(m3) + 0.5 * m3))

arms <- list(
  default = function() frm_sn(dd),
  raw_start = function() frm_sn(dd, start = raw_start),
  control_default = function() frm_sn(dd),
  arith = function() {
    s <- 0
    for (i in 1:2e5) s <- s + sqrt(i)
    s
  }
)
run <- function(f, reps) {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(reps)) f()
  proc.time()[["elapsed"]] - t0
}
# warmup and calibration: one fixed reps for every arm, from the
# slowest arm, so no arm's block is shorter than the next one's
for (nm in names(arms)) arms[[nm]]()
per <- vapply(arms, function(f) run(f, 20L) / 20, 0)
reps <- as.integer(ceiling(1.3 / max(per)))
cat("reps per block:", reps, " (slowest arm", format(max(per), digits = 3),
    "s per call)\n")

rounds <- 9L
nm_all <- names(arms)
m <- matrix(NA_real_, rounds, length(arms), dimnames = list(NULL, nm_all))
for (r in seq_len(rounds)) {
  ord <- nm_all[((seq_along(nm_all) + r - 2L) %% length(nm_all)) + 1L]
  for (nm in ord) m[r, nm] <- run(arms[[nm]], reps) / reps
}
best <- apply(m, 2, min)
cat("\nseconds per call, minimum of", rounds, "rounds\n")
for (nm in names(best)) cat(sprintf("  %-16s %10.5f\n", nm, best[[nm]]))
cat(sprintf("\nCONTROL  control_default / default  %.4f  (must be near 1)\n",
            best[["control_default"]] / best[["default"]]))
cat(sprintf("fallback raw_start / default          %.4f\n",
            best[["raw_start"]] / best[["default"]]))
cat(sprintf("arith control, cross-build load gauge %.5f s\n",
            best[["arith"]]))
esc <- frm_sn(dd, start = raw_start)$opt[["stationary_escape"]]
cat("raw_start arm escape:",
    if (is.null(esc)) "did not fire" else
      paste("fired,", esc[["starts"]], "restarts, gain",
            format(esc[["gain"]], digits = 6)), "\n")
