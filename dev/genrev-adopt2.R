.libPaths(c("C:/Users/adf44/source/r/genrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
f <- get("frm_adopt_generics", envir = asNamespace("frmtmb"))
block <- function(reps) { t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(reps)) f("frmtmb"); proc.time()[["elapsed"]] - t0 }
reps <- 64L
repeat { el <- block(reps); if (el >= 1.2 || reps > 2e6) break; reps <- reps * 4L }
best <- Inf; for (r in 1:5) best <- min(best, block(reps))
cat(sprintf("min-of-5 blocks: %d reps %.3f s -> %.1f us each\n",
            reps, best, 1e6 * best / reps))
cat(sprintf("single cold block of 1024: %.3f s -> %.1f us each\n",
            block(1024L), 1e6 * block(1024L) / 1024))
gc()
cat(sprintf("gc counts after: %s\n", paste(gc()[, "(Mb)"], collapse = ",")))
cat("GENREVDONE\n")
