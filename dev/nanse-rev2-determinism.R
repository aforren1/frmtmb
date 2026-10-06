# Reviewer, punch round 1 (m2): does the same frm() call give the same
# warnings? The fit-time Hessian budget reads wall-clock time
# (proc.time()[["elapsed"]] in fit_assembled() and fd_hessian()), so a
# fit near the budget may warn at fit time on one run and defer on
# another. Per model and replicate: did frm() itself warn, and was the
# check deferred (fit$cache$se_deferred set)?
#   Rscript dev/nanse-rev2-determinism.R <label> [reps]
args <- commandArgs(TRUE)
lab <- args[1]
reps <- if (length(args) > 1) as.integer(args[2]) else 10L
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
mk <- function(s) {
  set.seed(s)
  n <- 300
  x1 <- rnorm(n)
  g <- factor(rep(1:15, 20))
  data.frame(x1, g, y = 1 + 0.5 * x1 + rnorm(15, 0, 0.5)[g] + rnorm(n))
}
models <- list(
  ridge_re_s1 = list(d = mk(1), f = bf(y ~ a + c + dd, a ~ 0 + x1,
                                      c ~ 1 + (1 | g), dd ~ 1, nl = TRUE)),
  ridge_re_s11 = list(d = mk(11), f = bf(y ~ a + c + dd, a ~ 0 + x1,
                                        c ~ 1 + (1 | g), dd ~ 1, nl = TRUE))
)
for (nm in names(models)) {
  m <- models[[nm]]
  out <- character(reps)
  tf <- numeric(reps)
  for (r in seq_len(reps)) {
    w <- character()
    t0 <- proc.time()[["elapsed"]]
    f <- withCallingHandlers(frm(m$f, data = m$d),
      warning = function(x) {
        w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
      }, message = function(x) invokeRestart("muffleMessage"))
    tf[r] <- proc.time()[["elapsed"]] - t0
    se_w <- any(grepl("Standard errors are not available", w))
    out[r] <- if (se_w) "W" else if (!is.null(f$cache$se_deferred)) "D" else "-"
  }
  cat(sprintf("%s %-14s fit-time warning %2d, deferred %2d, neither %2d of %d | %s | elapsed %.3f to %.3f s\n",
              lab, nm, sum(out == "W"), sum(out == "D"), sum(out == "-"),
              reps, paste(out, collapse = ""), min(tf), max(tf)))
}
