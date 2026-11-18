# Short check of lane squant fix round 1, item 4 follow-up: the gp() fit
# of squant-short-compat.R ran for minutes. Time the same formula and data
# under gaussian() and quantile_elf(), with fewer rows, to tell a slow
# family from a slow gp().
# Run (worktree root):
#   Rscript dev/reviews/squant-short-gp.R > dev/reviews/squant-short-gp.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(20261117)
n <- 400
d <- data.frame(x = runif(n), x2 = runif(n),
                g = factor(sample(sprintf("g%02d", 1:25), n, TRUE)))
re0 <- rnorm(25, 0, 0.6)
re1 <- rnorm(25, 0, 0.5)
gi <- as.integer(d$g)
d$y <- 1 + 2 * d$x + re0[gi] + re1[gi] * d$x + (0.5 + d$x) * rnorm(n)
d$ys <- sin(2 * pi * d$x) + 0.5 * d$x2^2 + 0.4 * rnorm(n)
run <- function(lab, expr) {
  ws <- character()
  ms <- character()
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, substr(conditionMessage(w), 1, 150))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    ms <<- c(ms, substr(conditionMessage(m), 1, 150))
    invokeRestart("muffleMessage")
  }), error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(r, "error")) {
    cat(sprintf("%s: ERROR %s (%.1f s)\n", lab, conditionMessage(r), el))
    return(invisible(NULL))
  }
  cat(sprintf("%s: conv %s | logLik %.4f | warnings %d | messages %d | %.1f s\n",
              lab, r$opt$convergence, as.numeric(logLik(r)), length(ws),
              length(ms), el))
  if (length(ws)) cat("   W:", paste(unique(ws), collapse = "\n   W: "), "\n")
  if (length(ms)) cat("   M:", paste(unique(ms), collapse = "\n   M: "), "\n")
  invisible(r)
}
for (m in c(40, 80)) {
  dm <- d[seq_len(m), ]
  run(sprintf("gaussian gp(x) + x2, n %d", m),
      frm(ys ~ gp(x) + x2, family = gaussian(), data = dm))
  run(sprintf("quantile_elf gp(x) + x2, q 0.5, n %d", m),
      frm(bf(ys ~ gp(x) + x2, quantile = 0.5), family = quantile_elf(),
          data = dm))
  run(sprintf("quantile_elf(smoothing = 1) gp(x) + x2, q 0.5, n %d", m),
      frm(bf(ys ~ gp(x) + x2, quantile = 0.5),
          family = quantile_elf(smoothing = 1), data = dm))
}
run("gaussian gp(x) + x2, n 150",
    frm(ys ~ gp(x) + x2, family = gaussian(), data = d[1:150, ]))
cat("DONE\n")
